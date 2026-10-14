#!/usr/bin/env python3
"""
Python script to embed protein sequences from a FASTA file using ESM2 (esm2_t12_35M_UR50D).
 
For each sequence, the per-residue hidden states from the model's final
layer are mean-pooled (excluding special tokens) to produce a single
480-dimensional vector. The output is a tensor of shape (n_seqs, 480).
 
Usage:
    python esm2_embed.py input.fasta output.pt
    python esm2_embed.py input.fasta output.npy --format npy
"""
 
import argparse
import sys
 
import torch
from transformers import AutoTokenizer, AutoModel
 
 
MODEL_NAME = "facebook/esm2_t12_35M_UR50D"
HIDDEN_DIM = 480  # embedding dimension for this checkpoint
 
 
def read_fasta(path):
    """Yield (header, sequence) tuples from a FASTA file."""
    header = None
    seq_chunks = []
    with open(path, "r") as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)
 
 
@torch.no_grad()
def embed_sequences(headers, sequences, model, tokenizer, device, batch_size=8):
    """Return a (n_seqs, hidden_dim) tensor of mean-pooled embeddings."""
    all_embeddings = []
 
    for start in range(0, len(sequences), batch_size):
        batch_seqs = sequences[start:start + batch_size]
        batch_headers = headers[start:start + batch_size]
 
        encoded = tokenizer(
            batch_seqs,
            return_tensors="pt",
            padding=True,
            truncation=False,
        ).to(device)
 
        outputs = model(**encoded)
        hidden_states = outputs.last_hidden_state  # (B, L, hidden_dim)
        attention_mask = encoded["attention_mask"]  # (B, L)
 
        # ESM adds a BOS (<cls>) and EOS (<eos>) token. Build a mask that
        # zeroes out special tokens so they don't contribute to the mean.
        special_tokens_mask = tokenizer.get_special_tokens_mask(
            encoded["input_ids"][0].tolist(), already_has_special_tokens=True
        )
        # Recompute per-sequence since padding differs; simplest robust way:
        # mask = attention_mask AND NOT special-token position.
        input_ids = encoded["input_ids"]
        special_ids = {
            tokenizer.cls_token_id,
            tokenizer.eos_token_id,
            tokenizer.pad_token_id,
        }
        real_token_mask = torch.ones_like(attention_mask)
        for tok_id in special_ids:
            if tok_id is not None:
                real_token_mask &= (input_ids != tok_id)
        real_token_mask = real_token_mask * attention_mask  # (B, L)
 
        mask_expanded = real_token_mask.unsqueeze(-1).to(hidden_states.dtype)  # (B, L, 1)
        summed = (hidden_states * mask_expanded).sum(dim=1)  # (B, hidden_dim)
        counts = mask_expanded.sum(dim=1).clamp(min=1)  # (B, 1)
        mean_pooled = summed / counts  # (B, hidden_dim)
 
        all_embeddings.append(mean_pooled.cpu())
 
        for h in batch_headers:
            print(f"  embedded: {h}", file=sys.stderr)
 
    return torch.cat(all_embeddings, dim=0)
 
 
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_path", help="Path to input FASTA file")
    parser.add_argument(
        "output_path", help="Path to save the output tensor (.pt or .npy)"
    )
    parser.add_argument(
        "--format",
        choices=["pt", "npy", "csv"],
        default=None,
        help="Output format. Inferred from output_path extension if not set.",
    )
    parser.add_argument(
        "--batch-size", type=int, default=8, help="Batch size for inference"
    )
    args = parser.parse_args()
 
    fmt = args.format
    if fmt is None:
        if args.output_path.endswith(".npy"):
            fmt = "npy"
        elif args.output_path.endswith(".csv"):
            fmt = "csv"
        else:
            fmt = "pt"
 
    headers, sequences = zip(*read_fasta(args.fasta_path))
    headers, sequences = list(headers), list(sequences)
    print(f"Loaded {len(sequences)} sequences from {args.fasta_path}", file=sys.stderr)
 
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}", file=sys.stderr)
 
    print(f"Loading model {MODEL_NAME} ...", file=sys.stderr)
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = AutoModel.from_pretrained(MODEL_NAME).to(device)
    model.eval()
 
    embeddings = embed_sequences(
        headers, sequences, model, tokenizer, device, batch_size=args.batch_size
    )
    assert embeddings.shape == (len(sequences), HIDDEN_DIM), (
        f"Unexpected output shape {tuple(embeddings.shape)}, "
        f"expected ({len(sequences)}, {HIDDEN_DIM})"
    )
    print(f"Final embedding tensor shape: {tuple(embeddings.shape)}", file=sys.stderr)
 
    if fmt == "npy":
        import numpy as np
        np.save(args.output_path, embeddings.numpy())
    elif fmt == "csv":
        import csv
        arr = embeddings.numpy()
        with open(args.output_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["sequence_id"] + [f"dim_{i}" for i in range(arr.shape[1])])
            for header, row in zip(headers, arr):
                writer.writerow([header] + row.tolist())
    else:
        torch.save(embeddings, args.output_path)
 
    print(f"Saved embeddings to {args.output_path}", file=sys.stderr)
 
 
if __name__ == "__main__":
    main()