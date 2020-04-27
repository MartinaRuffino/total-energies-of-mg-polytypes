#% map kz=(q^2+.01^2-qx^2-qy^2)^.5*(qz>0?1:-1)
41  .5 .5 .5     0  0 0                L to Gamma
41   0  0  0     1  0 0                Gamma to X
21   1  0  0     1 .5 0                X to W
41   1 .5  0     0  0 0                W to Gamma
#41  0  0  0    .75 .75 0              Gamma to K
0    0 0 0  0 0 0

