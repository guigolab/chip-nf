#!/bin/bash
set -e
set -o pipefail

DB_FILE=chipseq-pipeline.db
VALIDATE_DIR=validate-ci

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
    ./nextflow run . $@
    ;;
  validate)
    shift
    f=$(getPath ${1-data/md5s})
    echo "Validating test results..." >&2
    [[ -s ${DB_FILE} ]] || false
    mkdir -p ${VALIDATE_DIR} && cd ${VALIDATE_DIR}
    cut -f 2 ../${DB_FILE} | xargs -I{} ln -fs {}
    md5sum -c ${f}
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
