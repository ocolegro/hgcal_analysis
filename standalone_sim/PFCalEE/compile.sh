#!/bin/bash
source g4env.sh && mkdir -p userlib/{lib,obj,bin} && cd userlib && make dictionary && make -j 5 && cd - && make -j 5


