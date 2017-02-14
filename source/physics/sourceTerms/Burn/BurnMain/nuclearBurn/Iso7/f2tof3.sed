s/use dBase, ONLY: ionmax/ /g
s/ionmax/NSPECIES/g
s/include 'eos_common.fh'/use Burn_dataEOS/g
s/include 'network_common.fh'/use Burn_data/g
s/rosen_gift/bn_rosenGift/g
s/rosen_ma28/bn_rosenMa28/g
s/stifbs_gift/bn_baderGift/g
s/stifbs3/bn_baderMa28/g
s/netint/bn_netIntegrate/g
s/initnet/bn_initNetwork/g
s/burner/bn_burner/g
s/azbar/bn_azbar/g
s/mazurek/bn_mazurek/g
s/screen4/bn_screen4/g
s/mcord/bn_mcord/g
s/sneutx/bn_sneutx/g
s/screen_iso7/bn_networkScreen/g
s/iso7tab/bn_networkTable/g
s/iso7rat/bn_networkRates/g
s/siso7/bn_bn_networkSparseJakob/g
s/biso7/bn_networkSparsePointers/g
s/diso7/bn_networkDenseJakob/g
s/iso7/bn_networkDerivs/g
s/use_table/bn_useBurnTable/g
s/ode_steper/bn_odeStepper/g


