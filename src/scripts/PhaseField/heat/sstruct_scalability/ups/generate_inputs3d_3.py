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
        <lower>[%(lower)f,%(lower)f,%(lower)f]</lower>
        <upper>[%(upper)f,%(upper)f,%(upper)f]</upper>
        <patches>[%(patches)d,%(patches)d,%(patches)d]</patches>
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
  <Solver type="hypre">
    <Parameters variable="u">
      <solver>%(solver)s</solver>
      <preconditioner>diagonal</preconditioner>
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

NAME = "sstruct_%(solver)s_3d_3_nlvl%(nlvl)1d_load%(load)04d"

for solver in ("fac","split","gmres","flexgmres","lgmres","bicgstab"):
	for nlvl in range(1,9):
		for pload in range(0,13,3):
			load = 2**pload;
			name = NAME % { "nlvl": nlvl, "load": load, "solver": solver };

			ups = open(name+".ups", "w")
			ups.write(HEAD % { "title": name })

			levels = """"""
			patches = 2**(1+pload/3)
			spacing = .4
			fullwidth = 16*patches*spacing*2**(nlvl-1);
			for lvl in range(1,nlvl+1):
				width = 16*patches*spacing;
				lower = (fullwidth-width)/2
				level = LEVEL % { "label": nlvl-lvl+1, "lower": lower, "upper": lower+width, "patches":patches, "spacing": spacing };
				levels = level + levels
				spacing = 2*spacing
				width = 3*width
			ups.write(levels)

			ups.write(FOOT % {'nlvl': nlvl, 'name': name, "solver": solver})
			ups.close()
