#!/bin/bash
export KMP_AFFINITY=verbose,granularity=core,compact
echo $KMP_AFFINITY
