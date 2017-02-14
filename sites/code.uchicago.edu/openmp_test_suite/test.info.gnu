<code>
  <UnitTest>
    <Eos>
      <Helmholtz>
        <3d>
          <UG>
            <twb>
              setupName: unitTest/Eos/Helmholtz
              setupOptions: -auto -test -3d +ug +noio threadWithinBlock=True -noc -makefile=gnu
              transfers: object/helm_table.dat object/SpeciesList.txt
              numProcs: 1
              parfiles: <pathToSimulations>/unitTest/Eos/Helmholtz/flash.par
            </twb>
          </UG>
          <nofbs>
            <twb>
              setupName: unitTest/Eos/Helmholtz
              setupOptions: -auto -test -3d +nofbs +noio threadWithinBlock=True -noc -makefile=gnu
              transfers: object/helm_table.dat object/SpeciesList.txt
              numProcs: 1
              parfiles: <pathToSimulations>/unitTest/Eos/Helmholtz/flash.par
            </twb>
          </nofbs>
        </3d>
      </Helmholtz>
    </Eos>
    <Gravity>
      <Poisson3>
        <3d>
          <UG>
            <twb>
              setupName: unitTest/Gravity/Poisson3
              setupOptions: -auto -test -3d +ug +cube32 +noio +newMpole threadWithinBlock=True -noc -makefile=gnu
              numProcs: 4
              parfiles: <pathToSimulations>/unitTest/Gravity/Poisson3/flash.par.ug
            </twb>
          </UG>
          <nofbs>
            <twb>
              setupName: unitTest/Gravity/Poisson3
              setupOptions: -auto -test -3d +nofbs +noio +newMpole threadWithinBlock=True -noc -makefile=gnu
              numProcs: 4
              parfiles: <pathToSimulations>/unitTest/Gravity/Poisson3/flash.par.ug
            </twb>
          </nofbs>
          <AMR>
            <tbl>
              setupName: unitTest/Gravity/Poisson3
              setupOptions: -auto -test -3d +pm4dev -maxblocks=600 +noio +newMpole threadBlockList=True -noc -makefile=gnu
              transfers: object/amr_runtime_parameters
              numProcs: 4
              parfiles: <pathToSimulations>/unitTest/Gravity/Poisson3/flash.par
            </tbl>
            <twb>
              setupName: unitTest/Gravity/Poisson3
              setupOptions: -auto -test -3d +pm4dev -maxblocks=600 +noio +newMpole -test threadWithinBlock=True -noc -makefile=gnu
              transfers: object/amr_runtime_parameters
              numProcs: 4
              parfiles: <pathToSimulations>/unitTest/Gravity/Poisson3/flash.par
            </twb>
          </AMR>
        </3d>
      </Poisson3>
    </Gravity>
  </UnitTest>
  <Comparison>
    <StirTurb>
      <AMR>
        <3d>
          <split>
            <tbl>
              setupName: StirTurb
              setupOptions: -auto -test -3d -unit=Particles +parallelio -noc threadBlockList=True -makefile=gnu
              numProcs: 4
              shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_StirTurb_AMR_3d_split/<runDir>/<chkMax>
              parfiles: <pathToSimulations>/StirTurb/coldstart_pm_3d.par
            </tbl>
            <twb>
              setupName: StirTurb
              setupOptions: -auto -test -3d -unit=Particles +parallelio -noc threadWithinBlock=True -makefile=gnu
              numProcs: 4
              shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_StirTurb_AMR_3d_split/<runDir>/<chkMax>
              parfiles: <pathToSimulations>/StirTurb/coldstart_pm_3d.par
            </twb>
          </split>
        </3d>
      </AMR>
      <nofbs>
        <3d>
          <split>
            <twb>
              setupName: StirTurb
              setupOptions: -auto -test -3d +nofbs -unit=Particles -noc -parfile=coldstart_nofbs_4p_3d.par threadWithinBlock=True -makefile=gnu
              numProcs: 4
              shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_StirTurb_nofbs_3d_split/<runDir>/<chkMax>
              parfiles: <pathToSimulations>/StirTurb/coldstart_nofbs_4p_3d.par
            </twb>
          </split>
        </3d>
      </nofbs>
    </StirTurb>
    <Sedov>
      <AMR>
        <2d>
          <split>
            <tbl>
              setupName: Sedov
              setupOptions: -auto -test +pm4dev +parallelio -noc threadBlockList=True -makefile=gnu
              transfers: object/amr_runtime_parameters
              numProcs: 4
              shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_Sedov_AMR_2d_split/<runDir>/<chkMax>
              parfiles: <pathToSimulations>/Sedov/coldstart_pm.par
            </tbl>
            <twb>
              setupName: Sedov
              setupOptions: -auto -test +pm4dev +parallelio -noc threadWithinBlock=True -makefile=gnu
              transfers: object/amr_runtime_parameters
              numProcs: 4
              shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_Sedov_AMR_2d_split/<runDir>/<chkMax>
              parfiles: <pathToSimulations>/Sedov/coldstart_pm.par
            </twb>
          </split>
          <unsplit>
            <tbl>
              setupName: Sedov
              setupOptions: -auto -test +uhd +pm4dev +parallelio -noc threadBlockList=True -makefile=gnu
              transfers: object/amr_runtime_parameters
              numProcs: 4
              shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_Sedov_AMR_2d_unsplit/<runDir>/<chkMax>
              parfiles: <pathToSimulations>/Sedov/coldstart_pm.par
            </tbl>
            <twb>
              setupName: Sedov
              setupOptions: -auto -test +uhd +pm4dev +parallelio -noc threadWithinBlock=True -makefile=gnu
              transfers: object/amr_runtime_parameters
              numProcs: 4
              shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_Sedov_AMR_2d_unsplit/<runDir>/<chkMax>
              parfiles: <pathToSimulations>/Sedov/coldstart_pm.par
            </twb>
          </unsplit>
        </2d>
      </AMR>
    </Sedov>
    <RTFlame>
      <nofbs>
        <2d>
          <split>
            <twb>
              setupName: RTFlame
              setupOptions: -auto -test +nofbs +parallelio -noc threadWithinBlock=True -makefile=gnu
              transfers: object/helm_table.dat
              numProcs: 4
              shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_RTFlame_nofbs_2d_split/<runDir>/<chkMax>
              parfiles: <pathToSimulations>/<setupName>/test_nofbs_2d.par
            </twb>
          </split>
        </2d>
      </nofbs>
      <AMR>
        <2d>
          <split>
            <tbl>
              setupName: RTFlame
              setupOptions: -auto -test +pm4dev -nxb=16 -nyb=16 +parallelio -noc threadBlockList=True -makefile=gnu
              transfers: object/amr_runtime_parameters object/helm_table.dat
              numProcs: 4
              shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_RTFlame_AMR_2d_split/<runDir>/<chkMax>
              parfiles: <pathToSimulations>/<setupName>/coldstart_pm_2d.par
            </tbl>
            <twb>
              setupName: RTFlame
              setupOptions: -auto -test +pm4dev -nxb=16 -nyb=16 +parallelio -noc threadWithinBlock=True -makefile=gnu
              transfers: object/amr_runtime_parameters object/helm_table.dat
              numProcs: 4
              shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_RTFlame_AMR_2d_split/<runDir>/<chkMax>
              parfiles: <pathToSimulations>/<setupName>/coldstart_pm_2d.par
            </twb>
          </split>
          <unsplit>
            <tbl>
              setupName: RTFlame
              setupOptions: -auto -test +uhd +pm4dev +parallelio -noc -nxb=16 -nyb=16 threadBlockList=True -makefile=gnu
              transfers: object/amr_runtime_parameters object/helm_table.dat
              numProcs: 4
              shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_RTFlame_AMR_2d_unsplit/<runDir>/<chkMax>
              parfiles: <pathToSimulations>/<setupName>/coldstart_pm_2d.par
            </tbl>
            <twb>
              setupName: RTFlame
              setupOptions: -auto -test +uhd +pm4dev +parallelio -noc -nxb=16 -nyb=16 threadWithinBlock=True -makefile=gnu
              transfers: object/amr_runtime_parameters object/helm_table.dat
              numProcs: 4
              shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_RTFlame_AMR_2d_unsplit/<runDir>/<chkMax>
              parfiles: <pathToSimulations>/<setupName>/coldstart_pm_2d.par
            </twb>
          </unsplit>
        </2d>
      </AMR>
    </RTFlame>
    <CurrentSheet>
      <AMR>
        <2d>
          <usm>
            <twb>
              setupName: magnetoHD/CurrentSheet
              setupOptions: -auto -test +pm4dev -2d +parallelio -noc +usm threadWithinBlock=True -makefile=gnu
              transfers: object/amr_runtime_parameters
              numProcs: 4
              shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_CurrentSheet_AMR_2d_usm/<runDir>/<chkMax>
              parfiles: <pathToSimulations>/<setupName>/coldstart_pm_2d.par
            </twb>
          </usm>
        </2d>
      </AMR>
    </CurrentSheet>
    <OrszagTang>
      <UG>
        <2d>
          <usm>
            <twb>
              setupName: magnetoHD/OrszagTang
              setupOptions: -auto -test +nofbs -2d +parallelio -noc +usm threadWithinBlock=True -makefile=gnu
              numProcs: 2
              shortPathToBenchmark: <siteDir>/gnu/2012-01-18/Comparison_OrszagTang_UG_2d_usm/<runDir>/<chkMax>
              parfiles: <pathToSimulations>/<setupName>/coldstart_nofbs_2p_2d.par
            </twb>
          </usm>
        </2d>
      </UG>
    </OrszagTang>
  </Comparison>
</code>
