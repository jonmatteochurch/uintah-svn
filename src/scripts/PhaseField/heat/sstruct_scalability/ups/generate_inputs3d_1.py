HEAD = """
<Uintah_specification>
  <Meta>
    <title>%(title)s</title>
  </Meta>
  <SimulationComponent type="phasefield"/>
  <PhaseField type="heat">
    <var>cc</var>
    <dim>3</dim>
    <delt>.2</delt>
    <alpha>1.</alpha>
    <refine_threshold>0.01</refine_threshold>
    <scheme>backward_euler</scheme>
  </PhaseField>
  <Time>
    <maxTime>20.2</maxTime>
    <initTime>0.</initTime>
    <delt_min>0.</delt_min>
    <delt_max>1.</delt_max>
    <timestep_multiplier>1.</timestep_multiplier>
  </Time>
  <Grid>"""

LEVEL = """
    <Level>
      <Box label="%(label)d">
        <lower>[0,0,0]</lower>
        <upper>[%(xupper)f,%(yupper)f,%(yupper)f]</upper>
        <patches>[%(xpatches)d,%(ypatches)d,%(ypatches)d]</patches>
      </Box>
      <spacing>[%(spacing)f,%(spacing)f,%(spacing)f]</spacing>
    </Level>"""

FOOT = """
    <BoundaryConditions>
      <Face side="x-">
        <BCType id="0" label="u" var="Neumann">
          <value>0.</value>
        </BCType>
      </Face>
      <Face side="y-">
        <BCType id="0" label="u" var="Neumann">
          <value>0.</value>
        </BCType>
      </Face>
      <Face side="z-">
        <BCType id="0" label="u" var="Neumann">
          <value>0.</value>
        </BCType>
      </Face>
      <Face side="x+">
        <BCType id="0" label="u" var="Dirichlet">
          <value>-0.65</value>
        </BCType>
      </Face>
      <Face side="y+">
        <BCType id="0" label="u" var="Dirichlet">
          <value>-0.65</value>
        </BCType>
      </Face>
      <Face side="z+">
        <BCType id="0" label="u" var="Dirichlet">
          <value>-0.65</value>
        </BCType>
      </Face>
    </BoundaryConditions>
  </Grid>
  <AMR>
    <Regridder type="Tiled">
      <adaptive>false</adaptive>
      <max_levels>%(nlvl)d</max_levels>
      <min_patch_size>[[16,16,16]]</min_patch_size>
      <min_boundary_cells>[1,1,1]</min_boundary_cells>
      <cell_refinement_ratio>[[2,2,2]]</cell_refinement_ratio>
      <cell_stability_dilation>[1,1,1]</cell_stability_dilation>
    </Regridder>
    <FineCoarseInterfaces>
      <FCIType id="1" label="u" var="FC1"/>
    </FineCoarseInterfaces>
  </AMR>
  <Solver type="hypre_sstruct" ndim="3" c2f="1">
    <Parameters variable="u">
      <solver>%(solver)s</solver>
      <setupFrequency>0</setupFrequency>
      <updateCoefFrequency>0</updateCoefFrequency>
      <solveFrequency>1</solveFrequency>
    </Parameters>
  </Solver>
  <DataArchiver>
    <filebase>%(name)s.uda</filebase>
    <outputTimestepInterval>0</outputTimestepInterval>
  </DataArchiver>
</Uintah_specification>"""

NAME = "sstruct_%(solver)s_3d_1_nlvl%(nlvl)1d_load%(load)04d"

for solver in ("fac","split","gmres","flexgmres","lgmres","bicgstab"):
	for nlvl in range(1,9):
		for pload in range(0,13,3):
			load = 2**pload;
			name = NAME % { "nlvl": nlvl, "load": load, "solver": solver };

			ups = open(name+".ups", "w")
			ups.write(HEAD % { "title": name })

			levels = """"""
			spacing = .4
			xpatches = 2**(1+pload/3)
			ypatches = 2**(nlvl+pload/3)
			yupper = 16*ypatches*spacing
			for lvl in range(1,nlvl+1):
				level = LEVEL % { "label": nlvl-lvl+1, "xupper": 16*xpatches*spacing, "yupper": yupper, "xpatches": xpatches, "ypatches": ypatches, "spacing": spacing };
				levels = level + levels
				spacing = spacing*2;
				ypatches = ypatches/2;
			ups.write(levels)

			ups.write(FOOT % {'nlvl': nlvl, 'name': name, "solver": solver})
			ups.close()
