#!/bin/bash

elements='ag al ar as au b ba be bi br c ca cd cl co cr cs cu f fe ga ge h he hf hg i in ir k kr li lu mg mn mo n na nb ne ni o os p pb pd po pt rb re rh rn ru s sb sc se si sn sr ta tc te ti tl v w xe y zn zr'

for i in $elements
do
	rm -f {log,mixm,moms,rst,save,atm,basp,wkp,site,ctrl,actrl,basp0,scflog*}.*$i
done
