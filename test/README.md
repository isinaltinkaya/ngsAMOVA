### Test


- test_mis, test_mis2

test_mis contains 14 sites: 9 missing for at least one individual and 5 that has data for both individuals.
test_mis2 contains only the 5 sites that has data for both individuals.

- test3.vcf

10 sites; (i1 missing sites: 2,6,8,10), (i2 missing sites: 8,10), (i3 missing sites: 4,8)
pairs; (i1,i2: 1,3,4,5,7,9) (i1,i3: 1,3,5,7,9) (i2,i3: 1,2,3,5,6,7,9)
onlyShared 1; (i1,i2: 1,3,5,7,9) (i1,i3: 1,3,5,7,9) (i2,i3: 1,3,5,7,9) 
onlyShared 0; (i1,i2: 1,3,4,5,7,9) (i1,i3: 1,3,5,7,9) (i2,i3: 1,2,3,5,6,7,9)
minInd 2 eq onlyShared 0; (i1,i2: 1,3,4,5,7,9) (i1,i3: 1,3,5,7,9) (i2,i3: 1,2,3,5,6,7,9)
minInd 3 eq onlyShared 1; (i1,i2: 1,3,5,7,9) (i1,i3: 1,3,5,7,9) (i2,i3: 1,3,5,7,9) 
