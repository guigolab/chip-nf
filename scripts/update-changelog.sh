#!/bin/bash
set -e
set -u
set -o pipefail

CHANGELOG_FILE="CHANGELOG.md"

make_changelog() {
    echo "$1"
    tail -n+2 $CHANGELOG_FILE
}

VER=$1
LOG=$(git log --format=format:"- %s" --grep 'ci skip' --grep 'skip ci' --invert-grep $(git describe --abbrev=0)...)

NEW_RELEASE=$(cat <<-CHANGELOG
# ChIP-nf Changelog

## Version ${VER}

$LOG
CHANGELOG)

make_changelog "$NEW_RELEASE" > .tmpcl && mv .tmpcl $CHANGELOG_FILE
