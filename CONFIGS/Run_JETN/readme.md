
# zonally periodic jet with internal waves configuration

The cpp key is JETN

The type of jet (1 to 4) is defined in `cppdefs.h`

Run time options are found in `croco.in`

# commands from start to finish:

```
./jobcomp
python chain_datarmor.py jet_cfg1_wp75_4km_0a2000j 20 03:00:00 4 jetn 0
# and release
```

Restart:
```
python chain_datarmor.py jet_cfg1_wp75_4km_nodecay_2000a3000j 10 03:00:00 4 jetn 1
cp /home2/pharos/othr/aponte/croco/jet_cfg1_wp75_4km_0a2000j/t20/jetn_rst*.nc /home1/scratch/aponte/jet_cfg1_wp75_4km_nodecay_2000a3000j/t0/
# and release
```



