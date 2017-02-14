<code>
  <Comparison>
    <StirTurb>
      <AMR>
        <3d>
          <split>
            setupName: StirTurb
            setupOptions: -auto -test -3d -unit=Particles +parallelio -noc -makefile=gnu
            numProcs: 4
            shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_StirTurb_AMR_3d_split/<runDir>/<chkMax>
            parfiles: <pathToSimulations>/StirTurb/coldstart_pm_3d.par
          </split>
        </3d>
      </AMR>
      <nofbs>
        <3d>
          <split>
            setupName: StirTurb
            setupOptions: -auto -test -3d +nofbs -unit=Particles -noc -parfile=coldstart_nofbs_4p_3d.par -makefile=gnu
            numProcs: 4
            shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_StirTurb_nofbs_3d_split/<runDir>/<chkMax>
            parfiles: <pathToSimulations>/StirTurb/coldstart_nofbs_4p_3d.par
          </split>
        </3d>
      </nofbs>
    </StirTurb>
    <Sedov>
      <AMR>
        <2d>
          <split>
            setupName: Sedov
            setupOptions: -auto -test +pm4dev +parallelio -noc -makefile=gnu
            transfers: object/amr_runtime_parameters
            numProcs: 4
            shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_Sedov_AMR_2d_split/<runDir>/<chkMax>
            parfiles: <pathToSimulations>/Sedov/coldstart_pm.par
          </split>
          <unsplit>
            setupName: Sedov
            setupOptions: -auto -test +uhd +pm4dev +parallelio -noc -makefile=gnu
            transfers: object/amr_runtime_parameters
            numProcs: 4
            shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_Sedov_AMR_2d_unsplit/<runDir>/<chkMax>
            parfiles: <pathToSimulations>/Sedov/coldstart_pm.par
          </unsplit>
        </2d>
      </AMR>
    </Sedov>
    <RTFlame>
      <nofbs>
        <2d>
          <split>
            setupName: RTFlame
            setupOptions: -auto -test +nofbs +parallelio -noc -makefile=gnu
            transfers: object/helm_table.dat
            numProcs: 4
            shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_RTFlame_nofbs_2d_split/<runDir>/<chkMax>
            parfiles: <pathToSimulations>/<setupName>/test_nofbs_2d.par
          </split>
        </2d>
      </nofbs>
      <AMR>
        <2d>
          <split>
            setupName: RTFlame
            setupOptions: -auto -test +pm4dev -nxb=16 -nyb=16 +parallelio -noc -makefile=gnu
            transfers: object/amr_runtime_parameters object/helm_table.dat
            numProcs: 4
            shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_RTFlame_AMR_2d_split/<runDir>/<chkMax>
            parfiles: <pathToSimulations>/<setupName>/coldstart_pm_2d.par
          </split>
          <unsplit>
	    setupName: RTFlame
            setupOptions: -auto -test +uhd +pm4dev +parallelio -noc -nxb=16 -nyb=16 -makefile=gnu
            transfers: object/amr_runtime_parameters object/helm_table.dat
            numProcs: 4
            shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_RTFlame_AMR_2d_unsplit/<runDir>/<chkMax>
            parfiles: <pathToSimulations>/<setupName>/coldstart_pm_2d.par
          </unsplit>
        </2d>
      </AMR>
    </RTFlame>
    <CurrentSheet>
      <AMR>
        <2d>
          <usm>
            setupName: magnetoHD/CurrentSheet
            setupOptions: -auto -test +pm4dev -2d +parallelio -noc +usm -makefile=gnu
            transfers: object/amr_runtime_parameters
            numProcs: 4
            shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_CurrentSheet_AMR_2d_usm/<runDir>/<chkMax>
            parfiles: <pathToSimulations>/<setupName>/coldstart_pm_2d.par
          </usm>
        </2d>
      </AMR>
    </CurrentSheet>
    <OrszagTang>
      <UG>
        <2d>
          <usm>
            setupName: magnetoHD/OrszagTang
            setupOptions: -auto -test +nofbs -2d +parallelio -noc +usm -makefile=gnu
            numProcs: 2
            shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_OrszagTang_UG_2d_usm/<runDir>/<chkMax>
            parfiles: <pathToSimulations>/<setupName>/coldstart_nofbs_2p_2d.par
          </usm>
        </2d>
      </UG>
    </OrszagTang>
  </Comparison>
</code>
