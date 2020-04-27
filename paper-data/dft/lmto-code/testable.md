The run path convention is program_system_other_stuff, the first _ is mandatory, the rest -- mere consistency.
To check times to split out the 'extra' tests:
  for i in `ninja test-show | grep testrun.sh | awk '{print $2}' | sed s:'./checks/'::g`; do echo $i; grep walltime checks/$i.log ; done
  for i in `ninja test-extra-show | grep testrun.sh | awk '{print $2}' | sed s:'./checks/'::g`; do echo $i; grep walltime checks/$i.log ; done
  for i in `ninja test-mpi-show | grep testrun.sh | awk '{print $2}' | sed s:'./checks/'::g`; do echo $i; grep walltime checks/$i.log ; done
  for i in `ninja test-mpi-extra-show | grep testrun.sh | awk '{print $2}' | sed s:'./checks/'::g`; do echo $i; grep walltime checks/$i.log ; done
With some exceptions, I'd put anything taking more than roughly 30s as 'extra'. Target 'test' then takes about 5min on an 8 core node.

| command                                                                    | run path           | n  | tags               |
|----------------------------------------------------------------------------|--------------------|----|--------------------|
| fp/test/test.fp --quiet --poszer          copt                             | fp_copt            |    |                    |
| fp/test/test.fp --quiet --poszer          te                               | fp_te              |    |                    |
| fp/test/test.fp --quiet --poszer          zrt                              | fp_zrt             |    |                    |
| fp/test/test.fp --quiet --poszer          co                               | fp_co              |    |                    |
| fp/test/test.fp --quiet --poszer          cr3si6                           | fp_cr3si6          |    |                    |
| fp/test/test.fp --quiet --poszer          fe                               | fp_fe              |    |                    |
| fp/test/test.fp --quiet --poszer          cu                               | fp_cu              |    |                    |
| fp/test/test.fp --quiet --poszer          au                               | fp_au              |    | extra              |
| fp/test/test.fp --quiet --poszer          gas                              | fp_gas             |    | extra              |
| fp/test/test.fp --quiet --poszer          srtio3                           | fp_srtio3          |    |                    |
| fp/test/test.fp --quiet --poszer          cdte 1                           | fp_cdte_1          |    |                    |
| fp/test/test.fp --quiet --poszer          cs 1                             | fp_cs_1            |    |                    |
| fp/test/test.fp --quiet --poszer          cdte 4                           | fp_cdte_4          |    |                    |
| fp/test/test.fp --quiet --poszer          felz                             | fp_felz            |    |                    |
| fp/test/test.fp --quiet --poszer          gasls                            | fp_gasls           |    |                    |
| fp/test/test.fp --quiet --poszer          zbgan                            | fp_zbgan           |    | extra              |
| fp/test/test.fp --quiet --poszer          gdn                              | fp_gdn             |    | extra              |
| fp/test/test.fp --quiet --poszer          eras                             | fp_eras            |    |                    |
| fp/test/test.fp --quiet --poszer          al                               | fp_al              |    |                    |
| fp/test/test.fp --quiet --poszer          c                                | fp_c               |    |                    |
| fp/test/test.fp --quiet --poszer          fept                             | fp_fept            |    | extra              |
| fp/test/test.fp --quiet --poszer          mgo                              | fp_mgo             |    | extra              |
| fp/test/test.fp --quiet --poszer          crn                              | fp_crn             |    | extra              |
|                                                                            |                    |    |                    |
| testing/test.lm --quiet --poszer          gas                              | lm_gas             |    |                    |
| testing/test.lm --quiet --poszer          tc                               | lm_tc              |    |                    |
| testing/test.lm --quiet --poszer          fe                               | lm_fe              |    |                    |
| testing/test.lm --quiet --poszer          mix                              | lm_mix             |    |                    |
| testing/test.lm --quiet --poszer          gan                              | lm_gan             |    |                    |
| testing/test.lm --quiet --poszer          cr3si6                           | lm_cr3si6          |    |                    |
| testing/test.lm --quiet --poszer          gd                               | lm_gd              |    |                    |
| testing/test.lm --quiet --poszer          er                               | lm_er              |    |                    |
| testing/test.lm --quiet --poszer          lap                              | lm_lap             |    |                    |
| testing/test.lm --quiet --poszer          si                               | lm_si              |    |                    |
| testing/test.lm --quiet --poszer          inp                              | lm_inp             |    |                    |
| testing/test.lm --quiet --poszer          feal                             | lm_feal            |    |                    |
| testing/test.lm --quiet --poszer          kfese                            | lm_kfese           |    | extra              |
|                                                                            |                    |    |                    |
| testing/test.blm --quiet --poszer --withlmfa          1                    | blm_1              |    |                    |
| testing/test.blm --quiet --poszer          2                               | blm_2              |    |                    |
| testing/test.blm --quiet --poszer          3                               | blm_3              |    |                    |
| testing/test.blm --quiet --poszer          4                               | blm_4              |    |                    |
| testing/test.blm --quiet --poszer          5                               | blm_5              |    |                    |
| testing/test.blm --quiet --poszer          6                               | blm_6              |    |                    |
| testing/test.blm --quiet --poszer          7                               | blm_7              |    |                    |
| testing/test.blm --quiet --poszer          8                               | blm_8              |    |                    |
| testing/test.blm --quiet --poszer          9                               | blm_9              |    |                    |
| testing/test.blm --quiet --poszer          10                              | blm_a              |    |                    |
| testing/test.blm --quiet --poszer          11                              | blm_b              |    |                    |
| testing/test.blm --quiet --poszer          12                              | blm_c              |    |                    |
|                                                                            |                    |    |                    |
| testing/test.lmscell --quiet --poszer          1                           | lmscell_1          |    |                    |
| testing/test.lmscell --quiet --poszer          2                           | lmscell_2          |    |                    |
| testing/test.lmscell --quiet --poszer          3                           | lmscell_3          |    |                    |
| testing/test.lmscell --quiet --poszer          4                           | lmscell_4          |    |                    |
| testing/test.lmscell --quiet --poszer          5                           | lmscell_5          |    |                    |
| testing/test.lmscell --quiet --poszer          6                           | lmscell_6          |    |                    |
| testing/test.lmscell --quiet --poszer          7                           | lmscell_7          |    |                    |
| testing/test.lmscell --quiet --poszer          8                           | lmscell_8          |    |                    |
|                                                                            |                    |    |                    |
| testing/test.ovlp --quiet --poszer          1                              | ovlp_1             |    |                    |
| testing/test.ovlp --quiet --poszer          2                              | ovlp_2             |    |                    |
| testing/test.ovlp --quiet --poszer          3                              | ovlp_3             |    |                    |
| testing/test.ovlp --quiet --poszer          4                              | ovlp_4             |    |                    |
| testing/test.ovlp --quiet --poszer          5                              | ovlp_5             |    |                    |
| testing/test.ovlp --quiet --poszer          6                              | ovlp_6             |    |                    |
| testing/test.ovlp --quiet --poszer          7                              | ovlp_7             |    |                    |
| testing/test.ovlp --quiet --poszer          8                              | ovlp_8             |    |                    |
| testing/test.ovlp --quiet --poszer          9                              | ovlp_9             |    |                    |
|                                                                            |                    |    |                    |
| testing/test.lmxbs --quiet --poszer          1                             | xbs_1              |    |                    |
|                                                                            |                    |    |                    |
| nc/test/test.nc --quiet --poszer          1                                | nc_1               |    |                    |
| nc/test/test.nc --quiet --poszer          2                                | nc_2               |    |                    |
| nc/test/test.nc --quiet --poszer          3                                | nc_3               |    |                    |
| nc/test/test.nc --quiet --poszer          4                                | nc_4               |    |                    |
| nc/test/test.nc --quiet --poszer          5                                | nc_5               |    |                    |
| nc/test/test.nc --quiet --poszer          6                                | nc_6               |    |                    |
| nc/test/test.nc --quiet --poszer          7                                | nc_7               |    |                    |
| nc/test/test.nc --quiet --poszer          8                                | nc_8               |    |                    |
|                                                                            |                    |    |                    |
| nc/test/test.so --quiet --poszer          1                                | nc_so1             |    |                    |
| nc/test/test.so --quiet --poszer          2                                | nc_so2             |    |                    |
|                                                                            |                    |    |                    |
| sx/test/test.sx --quiet --poszer          cdte 1                           | sx_cdte1           |    |                    |
| sx/test/test.sx --quiet --poszer          cdte 2                           | sx_cdte2           |    |                    |
| sx/test/test.sx --quiet --poszer          gan  1                           | sx_gan1            |    |                    |
|                                                                            |                    |    |                    |
| optics/test/test.optics --quiet --poszer          ogan 1                   | op_gan1            |    |                    |
| optics/test/test.optics --quiet --poszer          fe 1                     | op_fe1             |    |                    |
| optics/test/test.optics --quiet --poszer          mngas 1                  | op_mngas1          |    |                    |
| optics/test/test.optics --quiet --poszer          ogan 2                   | op_gan2            |    |                    |
| optics/test/test.optics --quiet --poszer          mngas 2                  | op_mngas2          |    |                    |
| optics/test/test.optics --quiet --poszer          ogan 3                   | op_gan3            |    |                    |
| optics/test/test.optics --quiet --poszer          fe 4                     | op_fe4             |    |                    |
| optics/test/test.optics --quiet --poszer          ogan 5                   | op_gan5            |    |                    |
| optics/test/test.optics --quiet --poszer          fe 5                     | op_fe5             |    |                    |
| optics/test/test.optics --quiet --poszer          ogan 6                   | op_gan6            |    |                    |
| optics/test/test.optics --quiet --poszer          fe 6                     | op_fe6             |    |                    |
| optics/test/test.optics --quiet --poszer          si 6                     | op_si6             |    |                    |
| optics/test/test.optics --quiet --poszer          fe 7                     | op_fe7             |    |                    |
| optics/test/test.optics --quiet --poszer          ogan 8                   | op_gan8            |    |                    |
| optics/test/test.optics --quiet --poszer          sic 9                    | op_sic9            |    |                    |
| optics/test/test.optics --quiet --poszer          ogan 10                  | op_gan10           |    |                    |
|                                                                            |                    |    |                    |
| gf/test/test.gf --quiet --poszer          co 1                             | gf_co1             |    |                    |
| gf/test/test.gf --quiet --poszer          co 2                             | gf_co2             |    |                    |
| gf/test/test.gf --quiet --poszer          co 4                             | gf_co4             |    |                    |
| gf/test/test.gf --quiet --poszer          co 5                             | gf_co5             |    |                    |
| gf/test/test.gf --quiet --poszer          co 7                             | gf_co7             |    |                    |
| gf/test/test.gf --quiet --poszer          gas 4                            | gf_gas4            |    |                    |
| gf/test/test.gf --quiet --poszer          fe 5                             | gf_fe5             |    |                    |
| gf/test/test.gf --quiet --poszer          fccfe 7                          | gf_fccfe7          |    |                    |
| gf/test/test.gf --quiet --poszer          mnpt 2                           | gf_mnpt2           |    |                    |
| gf/test/test.gf --quiet --poszer          mnpt 4                           | gf_mnpt4           |    |                    |
| gf/test/test.gf --quiet --poszer          mnpt 5                           | gf_mnpt5           |    |                    |
| gf/test/test.gf --quiet --poszer          mnpt 6                           | gf_mnpt6           |    |                    |
| gf/test/test.gf --quiet --poszer          mnn  1                           | gf_mnn1            |    |                    |
| gf/test/test.gf --quiet --poszer          mnn  2                           | gf_mnn2            |    |                    |
| gf/test/test.gf --quiet --poszer          mnn  4                           | gf_mnn4            |    |                    |
| gf/test/test.gf --quiet --poszer          mnn  5                           | gf_mnn5            |    |                    |
| gf/test/test.gf --quiet --poszer          mnn  6                           | gf_mnn6            |    |                    |
| gf/test/test.gf --quiet --poszer          ni   5                           | gf_ni5             |    |                    |
| gf/test/test.gf --quiet --poszer          cdte 8                           | gf_cdte8           |    |                    |
| gf/test/test.gf --quiet --poszer          eras 1                           | gf_eras1           |    |                    |
| gf/test/test.gf --quiet --poszer          eras 2                           | gf_eras2           |    |                    |
| gf/test/test.gf --quiet --poszer          eras 4                           | gf_eras4           |    |                    |
| gf/test/test.gf --quiet --poszer          eras 7                           | gf_eras7           |    |                    |
| gf/test/test.gf --quiet --poszer          fepd 9                           | gf_fepd9           |    |                    |
| gf/test/test.gf --quiet --poszer          fev  5                           | gf_fev             |    |                    |
| gf/test/test.gf --quiet --poszer          fev 10                           | gf_fev10           |    |                    |
| gf/test/test.frgf --quiet --poszer          ni 1                           | gf_frni1           |    |                    |
| gf/test/test.gf --quiet --poszer          femnpt  5                        | gf_femnpt5         |    | extra              |
| gf/test/test.gf --quiet --poszer          femnpt 10                        | gf_femnpt10        |    |                    |
| gf/test/test.gf --quiet --poszer          nife 6                           | gf_nife6           |    | extra              |
| gf/test/test.frgf --quiet --poszer          fept 2                         | gf_frfept2         |    | extra              |
|                                                                            |                    |    |                    |
| pgf/test/test.pgf --quiet --poszer          fe 1                           | pgf_fe1            |    |                    |
| pgf/test/test.pgf --quiet --poszer          fe 2                           | pgf_fe2            |    |                    |
| pgf/test/test.pgf --quiet --poszer          fe 3                           | pgf_fe3            |    |                    |
| pgf/test/test.pgf --quiet --poszer          fe 5                           | pgf_fe5            |    |                    |
| pgf/test/test.pgf --quiet --poszer          cr3si6 1                       | pgf_cr3si61        |    |                    |
| pgf/test/test.pgf --quiet --poszer          cr3si6 2                       | pgf_cr3si62        |    |                    |
| pgf/test/test.pgf --quiet --poszer          femgo 4                        | pgf_femgo4         |    |                    |
| pgf/test/test.pgf --quiet --poszer          --declead copt 5               | pgf_declead_copt5  |    |                    |
| pgf/test/test.pgf --quiet --poszer          --lu copt 5                    | pgf_lu_copt5       |    |                    |
| pgf/test/test.pgf --quiet --poszer          --lu=no copt 5                 | pgf_emb_copt5      |    |                    |
| pgf/test/test.pgf --quiet --poszer          --lu copt 7                    | pgf_lu_copt7       |    |                    |
| pgf/test/test.pgf --quiet --poszer          --lu=no copt 7                 | pgf_emb_copt7      |    |                    |
| pgf/test/test.pgf --quiet --poszer          --lu copt 8                    | pgf_lu_copt8       |    |                    |
| pgf/test/test.pgf --quiet --poszer          --lu=no copt 8                 | pgf_emb_copt8      |    |                    |
| pgf/test/test.pgf --quiet --poszer          --lu copt 9                    | pgf_lu_copt9       |    |                    |
| pgf/test/test.pgf --quiet --poszer          --lu=no copt 9                 | pgf_emb_copt9      |    |                    |
| pgf/test/test.pgf --quiet --poszer          --lu copt 10                   | pgf_lu_copt10      |    |                    |
| pgf/test/test.pgf --quiet --poszer          --lu=no copt 10                | pgf_emb_copt10     |    |                    |
| pgf/test/test.pgf --quiet --poszer          --lu copt 11                   | pgf_lu_copt11      |    |                    |
| pgf/test/test.pgf --quiet --poszer          --lu=no copt 11                | pgf_emb_copt11     |    |                    |
| pgf/test/test.pgf --quiet --poszer          --lu femgo 9                   | pgf_lu_femgo9      |    | extra              |
| pgf/test/test.pgf --quiet --poszer          --lu=no femgo 9                | pgf_emb_femgo9     |    | extra              |
| pgf/test/test.pgf --quiet --poszer          --lu femgo 11                  | pgf_lu_femgo11     |    | extra              |
| pgf/test/test.pgf --quiet --poszer          --lu=no femgo 11               | pgf_emb_femgo11    |    | extra              |
| pgf/test/test.pgf --quiet --poszer          co 2                           | pgf_co2            |    |                    |
| pgf/test/test.pgf --quiet --poszer          co 8                           | pgf_co8            |    | extra              |
| pgf/test/test.pgf --quiet --poszer          cuau 5                         | pgf_cuau5          |    |                    |
| pgf/test/test.pgf --quiet --poszer          free 8                         | pgf_free8          |    |                    |
| pgf/test/test.pgf --quiet --poszer          --lu coptco 7                  | pgf_lu_coptco7     |    | extra    disable   |
| pgf/test/test.pgf --quiet --poszer          --lu=no coptco 7               | pgf_emb_coptco7    |    | extra    disable   |
| pgf/test/test.pgf --quiet --poszer          --lu ava 10                    | pgf_lu_ava10       |    |                    |
| pgf/test/test.pgf --quiet --poszer          --lu=no ava 10                 | pgf_emb_ava10      |    |                    |
| pgf/test/test.pgf --quiet --poszer          aba 8                          | pgf_aba8           |    |                    |
| pgf/test/test.pgf --quiet --poszer          algas 9                        | pgf_algas9         |    | extra              |
|                                                                            |                    |    |                    |
| tb/test/test.tb --quiet --poszer          zrt                              | tb_zrt             |    |                    |
| tb/test/test.tb --quiet --poszer          tbso                             | tb_tbso            |    |          disable   |
| tb/test/test.tb --quiet --poszer          4h2o                             | tb_4h2o            |    |                    |
| tb/test/test.tb --quiet --poszer          fecr                             | tb_fecr            |    |                    |
| tb/test/test.tb --quiet --poszer          tbgan                            | tb_tbgan           |    |          disable   |
| tb/test/test.tb --quiet --poszer          tbovl                            | tb_tbovl           |    |          disable   |
| tb/test/test.tb --quiet --poszer          tbfit                            | tb_tbfit           |    |                    |
|                                                                            |                    |    |                    |
| mol/test/test.mol --quiet --poszer          h2o                            | mol_h2o            |    |                    |
| mol/test/test.mol --quiet --poszer          no2                            | mol_no2            |    |                    |
| mol/test/test.mol --quiet --poszer          h2op                           | mol_h2op           |    |                    |
|                                                                            |                    |    |                    |
| gwd/test/test.gwd --quiet --poszer          si                             | gwd_si             |    |                    |
| gwd/test/test.gwd --quiet --poszer          mno                            | gwd_mno            |    |                    |
| gwd/test/test.gwd --quiet --poszer          gas                            | gwd_gas            |    |                    |
| gwd/test/test.gwd --quiet --poszer          cu                             | gwd_cu             |    |                    |
| gwd/test/test.gwd --quiet --poszer          na                             | gwd_na             |    |                    |
|                                                                            |                    |    |                    |
| gwd/test/test.gwd --quiet --poszer          si2 2                          | gw_v0_si22         |    |                    |
| gwd/test/test.gwd --quiet --poszer          si2 4                          | gw_v0_si24         |    | extra              |
| gwd/test/test.gwd --quiet --poszer          --GWversion11 si2 2            | gw_v11_si22        |    |                    |
| gwd/test/test.gwd --quiet --poszer          --GWversion11 si2 4            | gw_v11_si24        |    | extra              |
| gwd/test/test.gwd --quiet --poszer          --GWversion12 si2 2            | gw_v12_si22        |    |                    |
| gwd/test/test.gwd --quiet --poszer          --GWversion12 si2 4            | gw_v12_si24        |    | extra              |
| gwd/test/test.gwd --quiet --poszer          fe 2                           | gw_fe2             |    | extra              |
| gwd/test/test.gwd --quiet --poszer          fe 4                           | gw_fe4             |    | extra              |
| gwd/test/test.gwd --quiet --poszer          coo2 2                         | gw_coo22           |    | extra              |
| gwd/test/test.gwd --quiet --poszer          coo2 4 7                       | gw_coo247          |    | extra              |
| gwd/test/test.gwd --quiet --poszer          cr3si6 7                       | gw_cr3si67         |    | extra              |
| gwd/test/test.gwd --quiet --poszer          six                            | gw_six             |    |                    |
| gwd/test/test.gwd --quiet --poszer          jell                           | gw_jell            |    |                    |
| gwd/test/test.nio6                                                         | gw_nio6            |    | extra              |
|                                                                            |                    |    |                    |
| dmft/test/test.dmft --quiet --poszer          lsco  1                      | dmft_lsco1         |    | extra              |
| dmft/test/test.dmft --quiet --poszer          lscoq 1                      | dmft_lscoq1        |    | extra              |
| dmft/test/test.dmft --quiet --poszer          lscoq 2                      | dmft_lscoq2        |    | extra              |
| dmft/test/test.dmft --quiet --poszer          ni 1                         | dmft_ni1           |    |                    |
| dmft/test/test.dmft --quiet --poszer          ni 2                         | dmft_ni2           |    | extra              |
| dmft/test/test.dmft --quiet --poszer          ni 3                         | dmft_ni3           |    | extra              |
| dmft/test/test.dmft --quiet --poszer          niso 1                       | dmft_niso1         |    |                    |
| dmft/test/test.dmft --quiet --poszer          niso 4                       | dmft_niso4         |    | extra              |
| dmft/test/test.dmft --quiet --poszer          nio 1                        | dmft_nio1          |    |                    |
|                                                                            |                    |    |                    |
| dmft/test/test.dmft --quiet --poszer --MPIK   ybco 1                       | dmft_ybco_mpik     | 6  |       mpik         |
| dmft/test/test_dmft_exchange_cix.py --mpi 2 --nk 2                         | dmft_exch_cix      | 2  | extra mpik         |
| dmft/test/test_dmft_gprt.py --mpi 2   --nk 2                               | dmft_gprt          | 2  |       mpik         |
| dmft/test/test_dmft_solver.py --mpi 2                                      | dmft_solver        | 2  |       mpik         |
| dmft/test/test_dmft_fullg.py --mpi 2                                       | dmft_fullg         | 2  |       mpik         |
|                                                                            |                    |    |                    |
| testing/test.scr --quiet 1                                                 | scr_1              |    |                    |
|                                                                            |                    |    |                    |
| testing/test.scr --quiet --poszer --MPIK                                   | scr_mpik           | 6  |       mpik         |
|                                                                            |                    |    |                    |
| testing/test.lm --quiet --poszer          --MPIK tc                        | lm_tc_mpik         | 6  |       mpik         |
| testing/test.lm --quiet --poszer          --MPIK gan                       | lm_gan_mpik        | 6  |       mpik         |
| testing/test.lm --quiet --poszer          --MPIK cr3si6                    | lm_cr3si6_mpik     | 8  |       mpik         |
| testing/test.lm --quiet --poszer          --MPIK gd                        | lm_gd_mpik         | 8  |       mpik         |
| testing/test.lm --quiet --poszer          --MPIK er                        | lm_er_mpik         | 8  |       mpik         |
| testing/test.lm --quiet --poszer          --MPIK si                        | lm_si_mpik         | 4  |       mpik         |
| testing/test.lm --quiet --poszer          --MPIK inp                       | lm_inp_mpik        | 5  |       mpik         |
| testing/test.lm --quiet --poszer          --MPIK feal                      | lm_feal_mpik       | 8  |       mpik         |
| testing/test.lm --quiet --poszer          --MPIK kfese                     | lm_kfese_mpik      | 6  |       mpik         |
|                                                                            |                    |    |                    |
| nc/test/test.nc --quiet --poszer          --MPIK 1                         | nc_1_mpik          | 8  |       mpik         |
| nc/test/test.nc --quiet --poszer          --MPIK 2                         | nc_2_mpik          | 8  |       mpik         |
| nc/test/test.nc --quiet --poszer          --MPIK 3                         | nc_3_mpik          | 8  |       mpik         |
| nc/test/test.nc --quiet --poszer          --MPIK 4                         | nc_4_mpik          | 8  |       mpik         |
| nc/test/test.nc --quiet --poszer          --MPIK 5                         | nc_5_mpik          | 8  |       mpik         |
| nc/test/test.nc --quiet --poszer          --MPIK 6                         | nc_6_mpik          | 8  |       mpik         |
| nc/test/test.nc --quiet --poszer          --MPIK 7                         | nc_7_mpik          | 8  |       mpik         |
| nc/test/test.nc --quiet --poszer          --MPIK 8                         | nc_8_mpik          | 8  |       mpik         |
|                                                                            |                    |    |                    |
| fp/test/test.fp --quiet --poszer          --MPIK co 1                      | fp_co1_mpik        | 4  |       mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK co 2                      | fp_co2_mpik        | 4  |       mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK copt 1                    | fp_copt_mpik       | 4  |       mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK te 1                      | fp_te_mpik         | 7  |       mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK zrt                       | fp_zrt_mpik        | 4  |       mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK cr3si6 1                  | fp_cr3si6_mpik     | 8  |       mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK cr3si6 2                  | fp_cr3si62_mpik    | 8  |       mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK fe 1                      | fp_fe1mpik         | 6  | extra mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK fe 2                      | fp_fe2mpik         | 6  | extra mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK cu 1                      | fp_cu_mpik         | 8  |       mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK gas 1                     | fp_gas_mpik        | 8  |       mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK srtio3 1                  | fp_srtio3_mpik     | 6  |       mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK cdte 1                    | fp_cdte1_mpik      | 4  |       mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK cdte 4                    | fp_cdte4_mpik      | 4  |       mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK felz 4                    | fp_felz_mpik       | 8  |       mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK zbgan 1                   | fp_zbgan_mpik      | 6  | extra mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK gdn 1                     | fp_gdn_mpik        | 4  |       mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK eras 1                    | fp_eras_mpik       | 7  |       mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK c                         | fp_c_mpik          | 4  |       mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK fept                      | fp_fept_mpik       | 6  |       mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK mgo 3                     | fp_mgo_mpik        | 7  |       mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK crn 2                     | fp_crn_mpik        | 8  |       mpik         |
|                                                                            |                    |    |                    |
| gf/test/test.gf --quiet --poszer          --MPIK gas                       | gf_gas_mpik        | 4  |       mpik         |
! gf/test/test.gf --quiet --poszer --MPIK co 1                               | gf_co1_mpik        | 4  |       mpik         |
| gf/test/test.gf --quiet --poszer --MPIK co 2                               | gf_co2_mpik        | 4  |       mpik         |
| gf/test/test.gf --quiet --poszer --MPIK co 4                               | gf_co4_mpik        | 4  |       mpik         |
| gf/test/test.gf --quiet --poszer --MPIK co 7                               | gf_co7_mpik        | 4  |       mpik         |
| gf/test/test.gf --quiet --poszer --MPIK gas 4                              | gf_gas4_mpik       | 4  |       mpik         |
| gf/test/test.gf --quiet --poszer --MPIK fccfe 7                            | gf_fccfe7_mpik     | 4  |       mpik         |
| gf/test/test.gf --quiet --poszer --MPIK mnpt 2                             | gf_mnpt2_mpik      | 4  |       mpik         |
| gf/test/test.gf --quiet --poszer --MPIK mnpt 4                             | gf_mnpt4_mpik      | 4  |       mpik         |
| gf/test/test.gf --quiet --poszer --MPIK mnn  1                             | gf_mnn1_mpik       | 4  |       mpik         |
| gf/test/test.gf --quiet --poszer --MPIK mnn  2                             | gf_mnn2_mpik       | 4  |       mpik         |
| gf/test/test.gf --quiet --poszer --MPIK mnn  4                             | gf_mnn4_mpik       | 4  |       mpik         |
| gf/test/test.gf --quiet --poszer --MPIK cdte 8                             | gf_cdte8_mpik      | 4  |       mpik         |
| gf/test/test.gf --quiet --poszer --MPIK eras 1                             | gf_eras1_mpik      | 4  |       mpik         |
| gf/test/test.gf --quiet --poszer --MPIK eras 2                             | gf_eras2_mpik      | 4  |       mpik         |
| gf/test/test.gf --quiet --poszer --MPIK eras 4                             | gf_eras4_mpik      | 4  |       mpik         |
| gf/test/test.gf --quiet --poszer --MPIK eras 7                             | gf_eras7_mpik      | 4  |       mpik         |
| gf/test/test.gf --quiet --poszer --MPIK fepd 9                             | gf_fepd9_mpik      | 4  |       mpik         |
| gf/test/test.gf --quiet --poszer --MPIK fev 10                             | gf_fev10_mpik      | 4  |       mpik         |
| gf/test/test.gf --quiet --poszer --MPIK femnpt 10                          | gf_femnpt10_mpik   | 4  |       mpik         |
|                                                                            |                    |    |                    |
| pgf/test/test.pgf --quiet --poszer          --MPIK copt 5                  | pgf_copt5_mpik     | 4  |       mpik         |
| pgf/test/test.pgf --quiet --poszer          --MPIK copt 7                  | pgf_copt7_mpik     | 4  |       mpik         |
| pgf/test/test.pgf --quiet --poszer          --MPIK copt 9                  | pgf_copt9_mpik     | 4  |       mpik         |
| pgf/test/test.pgf --quiet --poszer          --MPIK copt 10                 | pgf_copt10_mpik    | 4  |       mpik         |
| pgf/test/test.pgf --quiet --poszer          --MPIK copt 11                 | pgf_copt11_mpik    | 4  |       mpik         |
| pgf/test/test.pgf --quiet --poszer          --lu --MPIK femgo 9            | pgf_lu_femgo9_mpik | 4  | extra mpik         |
| pgf/test/test.pgf --quiet --poszer          --lu=no --MPIK femgo           | pgf_em_femgo9_mpik | 4  | extra mpik         |
| pgf/test/test.pgf --quiet --poszer          --MPIK femgo 11                | pgf_femgo11_mpik   | 4  | extra mpik         |
| pgf/test/test.pgf --quiet --poszer  --lu=no --MPIK coptco 7                | pgf_em_coptco7_mpik| 4  | extra mpik         |
| pgf/test/test.pgf --quiet --poszer     --lu --MPIK coptco 7                | pgf_coptco7_mpik   | 4  | extra mpik         |
|                                                                            |                    |    |                    |
| gwd/test/test.gwd --quiet --poszer          --mpi=4 si                     | gwd_si_mpik        | 4  |       mpik         |
| gwd/test/test.gwd --quiet --poszer          --mpi=4 mno                    | gwd_mno_mpik       | 4  |       mpik         |
| gwd/test/test.gwd --quiet --poszer          --mpi=4 gas                    | gwd_gas_mpik       | 4  |       mpik         |
| gwd/test/test.gwd --quiet --poszer          --mpi=4 cu                     | gwd_cu_mpik        | 4  |       mpik         |
| gwd/test/test.gwd --quiet --poszer          --mpi=4 na                     | gwd_na_mpik        | 4  |       mpik         |
|                                                                            |                    |    |                    |
| gwd/test/test.gwd --quiet --poszer --mpi=4,4 zbmnas 6                      | gw_zbmnas6_mpik    | 4  |       mpik         |
| gwd/test/test.gwd --quiet --poszer --mpi=1,4 lif 4 6                       | gw_lif_mpik        | 4  |       mpik         |
| gwd/test/test.gwd --quiet --poszer --mpi=4,4 six                           | gw_six_mpik        | 4  |       mpik         |
| gwd/test/test.gwd --quiet --poszer --mpi=4,4 jell                          | gw_jell_mpik       | 4  |       mpik         |
|                                                                            |                    |    |                    |
| gwd/test/test.gwd --quiet --poszer          --mpi=4,4 si2                  | gw_si2_mpik        | 4  | extra mpik         |
| gwd/test/test.gwd --quiet --poszer          --mpi=4,4 fe                   | gw_fe_mpik         | 4  | extra mpik         |
| gwd/test/test.gwd --quiet --poszer          --mpi=4,4 coo2                 | gw_coo2_mpik       | 4  | extra mpik         |
| gwd/test/test.gwd --quiet --poszer          --mpi=4,4 nio                  | gw_nio_mpik        | 4  | extra mpik         |
| gwd/test/test.gwd --quiet --poszer          --mpi=8,8 cr3si6               | gw_cr3si6_mpik     | 8  | extra mpik         |
|                                                                            |                    |    |                    |
| gw/test/test.gw --quiet --poszer          --mpi=8,8 cr3si6                 | gwx_cr3si6_mpik     | 8  | extra mpik        |

| testing/test.lm --quiet --poszer          --libxc kfese                    | lm_lxckfese        |    | extra              |
|                                                                            |                    |    |                    |
| fp/test/test.tio2                                                          | fp_tio2            |    | extra              |
| fp/test/test.fp --quiet --noplot --poszer na                               | fp_na              |    |                    |
| fp/test/test.fp --quiet --noplot --poszer coptso                           | fp_coptso          |    | extra      disable |
| fp/test/test.fp --quiet --noplot --poszer er                               | fp_er              |    | extra              |
| fp/test/test.fp --quiet --noplot --poszer kfese                            | fp_kfese           |    | extra              |
| fp/test/test.fp --quiet --noplot --poszer ni                               | fp_ni              |    |                    |
| fp/test/test.fp --quiet --noplot --poszer bzt                              | fp_bzt             |    | extra              |
|                                                                            |                    |    |                    |
| fp/test/test.fp --libxc --quiet --poszer te                                | fp_te_libxc        |    |                    |
| fp/test/test.fp --libxc --quiet --poszer fe 1                              | fp_fe_libxc        |    | extra              |
|                                                                            |                    |    |                    |
| gf/test/test.gf --quiet --poszer          fept 5                           | gf_fept-5          |    |                    |
| gf/test/test.gf --quiet --poszer          fept 7                           | gf_fept-7          |    |                    |
| gf/test/test.gf --quiet --poszer          fept2 7                          | gf_fept2-7         |    |                    |
| gf/test/test.gf --quiet --poszer          --codecheck fe2b                 | gf_fe2b            |    | extra      disable |
| gf/test/test.gf --quiet --poszer          fevg                             | gf_fevg            |    |                    |
| gf/test/test.frgfmst                                                       | gf_frgfmst         |    | extra              |
| gf/test/test.ncgfso                                                        | gf_ncgfso          |    |                    |
| gf/test/test.ncgfmst                                                       | gf_ncgfmst         |    |                    |
| gf/test/test.ncgfmst2                                                      | gf_ncgfmst2        |    |                    |
| gf/test/test.srgf                                                          | gf_srgf            |    | extra              |
| gf/test/test.frgfm1                                                        | gf_frgfm1          |    | extra              |
| gf/test/test.frgfm10                                                       | gf_frgfm10         |    |                    |
| gf/test/test.gf --quiet --gamma=4 femnpt 10                                | gf_femnpt-10       |    |                    |
|                                                                            |                    |    |                    |
| pgf/test/test.pgf --quiet --poszer --lu=no copt                            | pgf_copt_lu        |    | extra              |
| pgf/test/test.pgf --quiet --poszer --lu    fe 1                            | pgf_fe-lu1         |    |                    |
| pgf/test/test.pgf --quiet --poszer --lu    fe 5                            | pgf_fe-lu5         |    |                    |
| pgf/test/test.pgf --quiet --poszer --lu    aba                             | pgf_aba_lu         |    |                    |
| pgf/test/test.pgf --quiet --poszer --lu    ava                             | pgf_ava_lu         |    | extra              |
| pgf/test/test.pgf --quiet --poszer --MPIK --mergepl coptco 7               | pgf_coptcom_mpik   | 4  | extra mpik         |
|                                                                            |                    |    |                    |
| testing/test.lm --quiet --poszer          --MPIK --libxc kfese             | lm_lxckfese_mpik   | 6  | extra mpik         |
|                                                                            |                    |    |                    |
| dmft/test/test.dmft --quiet --poszer          --MPIK lscoq 1               | dmft_lscoq1_mpik   | 6  | extra mpik         |
| dmft/test/test.dmft --quiet --poszer          --MPIK lscoq 2               | dmft_lscoq2_mpik   | 6  | extra mpik         |
| dmft/test/test.dmft --quiet --poszer          --MPIK ni 1                  | dmft_ni1_mpik      | 6  | extra mpik         |
| dmft/test/test.dmft --quiet --poszer          --MPIK ni 2                  | dmft_ni2_mpik      | 6  | extra mpik         |
| dmft/test/test.dmft --quiet --poszer          --MPIK niso 1                | dmft_niso1_mpik    | 6  | extra mpik         |
| dmft/test/test.dmft --quiet --poszer          --MPIK nio 1                 | dmft_nio1_mpik     | 6  | extra mpik         |
| dmft/test/test.dmft --quiet --poszer          --MPIK niso 4                | dmft_niso4_mpik    | 6  | extra mpik         |
| dmft/test/test.dmft --quiet --poszer          --MPIK ni 3                  | dmft_ni3_mpik      | 6  | extra mpik         |
|                                                                            |                    |    |                    |
| gf/test/test.gf --quiet --MPIK fept                                        | gf_fept_mpik       | 4  |       mpik         |
| gf/test/test.gf --quiet --MPIK --codecheck fe2b                            | gf_fe2bcch_mpik    | 4  | extra mpik         |
|                                                                            |                    |    |                    |
| fp/test/test.tio2 --MPIK                                                   | fp_tio2_mpik       | 8  |       mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK coptso                    | fp_coptso_mpik     | 4  | extra mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK kfese                     | fp_kfese_mpik      | 6  | extra mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK bzt                       | fp_bzt_mpik        | 5  | extra mpik         |
| fp/test/test.fp --quiet --poszer          --MPIK er                        | fp_er_mpik         | 8  | extra mpik         |
|                                                                            |                    |    |                    |
| gwd/test/test.gwd --quiet --timereversaloff --eibzmodeoff coo2 2           | gw_tr_coo22        | 1  | extra              |
| gwd/test/test.gwd --quiet --timereversaloff --eibzmodeoff coo2 4           | gw_tr_coo24        | 1  | extra              |
|                                                                            |                    |    |                    |
| gwd/test/test.gwd --quiet --mpi=4,4 fe 2                                   | gw_fe2_mpik        | 4  | extra mpik         |
| gwd/test/test.gwd --quiet --mpi=4,4 fe 4 5                                 | gw_fe45_mpik       | 4  | extra mpik         |
| gwd/test/test.gwd --quiet --mpi=4,4 --so2 fe 4                             | gw_fe4_so2_mpik    | 4  | extra mpik         |
| gwd/test/test.gwd --quiet --mpi=4,4 --so3 fe 4                             | gw_fe4_so3_mpik    | 4  | extra mpik         |
| gwd/test/test.gwd --quiet --mpi=4,4 copt 2                                 | gw_copt2_mpik      | 4  | extra mpik         |
| gwd/test/test.gwd --quiet --mpi=4,4 copt 4                                 | gw_copt4_mpik      | 4  | extra mpik         |
| gwd/test/test.gwd --quiet --mpi=4,4 --eibzmodeoff copt 2 4                 | gw_copt2nosym_mpik | 4  | extra mpik         |
| gwd/test/test.gwd --quiet --mpi=4,4 --eibzmodeoff copt 4                   | gw_copt4nosym_mpik | 4  | extra mpik         |
| gwd/test/test.gwd --quiet --mpi=4,4 --timereversaloff copt 4               | gw_copt4tr_mpik    | 4  | extra mpik         |
| gwd/test/test.gwd --quiet --mpi=4,4 --timereversaloff --eibzmodeoff copt 4 | gw_copt4trns_mpik  | 4  | extra mpik         |
| gwd/test/test.gwd --quiet --mpi=4,4 nio 4                                  | gw_nio4_mpik       | 4  | extra mpik         |
| gwd/test/test.gwd --quiet --mpi=4,4 nio 6                                  | gw_nio6_mpik       | 4  | extra mpik         |
| gwd/test/test.gwd --quiet --eibzmodeoff --offbz --mpi=4,4 nio 6            | gw_niooffbz_mpik   | 4  | extra mpik         |
