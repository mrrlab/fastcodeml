#!/bin/bash
ls *.cpp *.h | grep -Ev "^(blas[.]h|lapack[.]h|nlopt[.]h)$" | xargs clang-format -i
