#!/bin/bash
set -e
set -o pipefail

DB_FILE=chipseq-pipeline.db
VALIDATE_DIR=validate-ci

export NXF_ANSI_LOG=false
export NXF_DOCKER=${NXF_DOCKER-'with'}
getPath() {
    if [ -n "$(type -a realpath)" ]; then
        realpath $@
    else
        readlink -f $@
    fi
}

case "$1" in
  run)
    shift
    echo "Running test pipeline..." >&2
    ./nextflow run . -${NXF_DOCKER}-docker $@
    ;;
  validate)
    shift
    md5s=$(getPath ${1-data/md5s})
    echo "Validating test results..." >&2
    [[ -s ${DB_FILE} ]] || false
    mkdir -p ${VALIDATE_DIR} && cd ${VALIDATE_DIR}
    cut -f 2 ../${DB_FILE} | while read f
    do
        name=$(basename $f)
        ext=${name##*.}
        case $ext in
            bam)
                samtools view ${f} > ${name%.*}.sam
                ;;
            *)
                ln -fs ${f} .
        esac
    done
    md5sum -c ${md5s}
    ;;
  cleanup)
    echo "Cleaning up test results..." >&2
    find ${VALIDATE_DIR} -type l -exec rm {} \+
    find ${VALIDATE_DIR} -type d -empty -exec rmdir {} \+
    ;;
  *)
    echo "Usage: ci {run|validate|cleanup}" >&2
    exit 1
esac
