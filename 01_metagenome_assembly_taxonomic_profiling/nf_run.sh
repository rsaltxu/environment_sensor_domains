#!/bin/bash

nextflow run nf-core/mag -profile conda -params-file nf-params.json -c base.config -resume