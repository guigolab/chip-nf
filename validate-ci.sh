#!/bin/bash
set -e
set -o pipefail

DB_FILE=chipseq-pipeline.db

[[ -s ${DB_FILE} ]] || false

DIR=validate-ci
mkdir -p ${DIR} && cd ${DIR}
cut -f 2 ../${DB_FILE} | xargs -I{} ln -fs {}
md5sum -c ../data/md5s
cd - && rm -rf ${DIR}