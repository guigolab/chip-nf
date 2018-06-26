# ChIP-nf Changelog

## Version 0.2.2

- Only specify -p option in model when cpus > 1 - fix #6
- Move to Travis CI

## Version 0.2.1

- Fix duplicate input file names in zerone - close #5
- Remove debugging code making release 0.2.0 invalid - #4

## Version 0.2.0

- Update groupTuple for merging and fix input sorting for Zerone - resolves #3
- Add paramter to set MACS2 temp dir and set default to process folder
- Use global fragment length from cli option when not shifting - close #2
- Fix issue with samples using the same control
- Add more test data
- Update pipeline with Zerone process and change config and readme

## Version 0.1.0

First version
