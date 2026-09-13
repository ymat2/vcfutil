#!/bin/bash

vcfutil -h
vcfutil daf -h

vcfutil daf --vcf tests/data/sample.vcf --out sample --window_size 100 --window_step 50 --pop1 tests/data/pop1.txt --pop2 tests/data/pop2.txt
vcfutil daf --vcf tests/data/sample.vcf --out sample --window_size 100 --pop1 tests/data/pop1.txt --pop2 tests/data/pop2.txt
vcfutil daf --vcf tests/data/sample.vcf --out sample --window_size 100 --window_step 50 --site 0.5 --pop1 tests/data/pop1.txt --pop2 tests/data/pop2.txt

vcfutil daf --vcf tests/data/sample.vcf --site 0.5 --pop1 tests/data/pop1.txt --pop2 tests/data/pop2.txt
vcfutil daf --vcf tests/data/sample.vcf --site2 0.6 --pop1 tests/data/pop1.txt --pop2 tests/data/pop2.txt
vcfutil daf --vcf tests/data/sample.vcf --site 0 --pop1 tests/data/pop1.txt --pop2 tests/data/pop2.txt
vcfutil daf --vcf tests/data/sample.vcf --site2 0 --pop1 tests/data/pop1.txt --pop2 tests/data/pop2.txt
vcfutil daf --vcf tests/data/sample.vcf --site 0 --site2 0.3 --pop1 tests/data/pop1.txt --pop2 tests/data/pop2.txt
