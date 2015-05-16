#!/bin/bash
#./fast-full --debug 2 --number-of-threads 1 --branch-lengths-fixed --branch 2 --no-pre-stop --output EMGT00050000008747.Drosophila.002.fcml EMGT00050000008747.Drosophila.002.nwk EMGT00050000008747.Drosophila.002.phy
./fast-aggr --debug 2 --number-of-threads 1 --branch-lengths-fixed --branch 2 --no-pre-stop --output EMGT00050000008747.Drosophila.002.fcml EMGT00050000008747.Drosophila.002.nwk EMGT00050000008747.Drosophila.002.phy
~/work/dndstools/fcml2tab.py --long --output EMGT00050000008747.Drosophila.002.tab EMGT00050000008747.Drosophila.002.fcml
./fast_cont.py --exe ./fast-full EMGT00050000008747.Drosophila.002.tab EMGT00050000008747.Drosophila.002.nwk EMGT00050000008747.Drosophila.002.phy EMGT00050000008747.Drosophila.002.fcont -- --branch-lengths-fixed --number-of-threads 1 --debug 2
