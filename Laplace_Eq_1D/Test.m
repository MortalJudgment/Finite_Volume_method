% Simulation parameters.
 nlev   = 1:3;        % Grid levels.
 nsolv  = 1:4;        % Solvers 1: FEATool, 2: FEniCS, 3: OpenFOAM, 4: SU2
 iplot  = 1;          % Plot and visualize results.
 isave  = true;       % Save benchmark results to file.
 iimg   = false;      % Save images.
 ldir   = fullfile(tempdir(),'ns3_bench');   % Benchmark directory.
 lfile  = 'ns3_bench.log';   % Logfile name.

% Test case settings.
 rho       = 1;           % Density.
 miu       = 0.001;       % Molecular/dynamic viscosity.
 umax      = 0.3;         % Maximum magnitude of inlet velocity.
 umean     = 2/3*umax;    % Mean inlet velocity.
% Geometry.
 h         = 0.41;        % Height of rectangular domain.
 l         = 2.2;         % Length of rectangular domain.
 xc        = 0.2;         % x-coordinate of cylinder center.
 yc        = 0.2;         % y-coordinate of cylinder center.
 diam      = 0.1;         % Diameter of cylinder.
% FEM discretization parameters.
 sf_u      = 'sflag2';    % FEM shape function type for velocity.
 sf_p      = 'sflag1';    % FEM shape function type for pressure.


% Run benchmarks.
 rmkdir(ldir);
 solvers = {'FEATool','FEniCS','OpenFOAM','SU2'};
 data = {};
 for isolv=nsolv
   data_i = [];
   for ilev=nlev
     clear fea

% Geometry definition.
     gobj1 = gobj_rectangle( 0, 2.2, 0, 0.41, 'R1' );
     gobj2 = gobj_circle( [0.2 0.2], 0.05, 'C1' );
     fea.geom.objects = { gobj1 gobj2 };
     fea = geom_apply_formula( fea, 'R1-C1' );
     fea.sdim = { 'x', 'y' };

% Grid generation.
     if( ilev>=1 )

% Structured quadrilateral benchmark grid.
       ns = 8*2^(ilev-1);
       r  = [0.05 0.06 0.08 0.11 0.15];
       x  = [0.41 0.5 0.7 1 1.4 1.8 2.2];
       for i=2:ilev
         r = sort( [ r, (r(1:end-1)+r(2:end))/2 ] );
         x = sort( [ x, (x(1:end-1)+x(2:end))/2 ] );
       end

       grid1 = ringgrid( r, 4*ns, [], [], [0.2;0.2] );
       grid2 = holegrid( ns, 2^(ilev-1), [0 0.41;0 0.41], 0.15, [0.2;0.2] );
       grid2 = gridmerge( grid1, 5:8, grid2, 1:4 );
       grid3 = rectgrid( x, ns, [0.41 2.2;0 0.41] );
       fea.grid = gridmerge( grid3, 4, grid2, 6 );
       fea.grid.s(:) = 1;

       if( isolv==2 && size(fea.grid.c,1)==4 )
         fea.grid = quad2tri( fea.grid );
       end
     else
       fea.grid = gridgen( fea, 'hmax', abs(ilev) );
     end

% Boundary specifications.
     DTOL      = sqrt(eps)*1e3;
     i_inflow  = findbdr( fea, ['x<=',num2str(DTOL)] );     % Inflow boundary number.
     i_outflow = findbdr( fea, ['x>=',num2str(l-DTOL)] );   % Outflow boundary number.
     s_inflow  = ['4*',num2str(umax),'*(y*(',num2str(h),'-y))/',num2str(h),'^2'];   % Definition of inflow profile.
     i_cyl     = findbdr( fea, ['sqrt((x-',num2str(xc),').^2+(y-',num2str(yc),').^2)<=(',num2str(diam/2+DTOL),')'] );    % Cylinder boundary number.

% Problem definition.
     fea = addphys(fea,@navierstokes);
     fea.phys.ns.eqn.coef{1,end} = { rho };
     fea.phys.ns.eqn.coef{2,end} = { miu };
     fea.phys.ns.sfun            = { sf_u, sf_u, sf_p };

% Boundary conditions.
     fea.phys.ns.bdr.sel(i_inflow)  = 2;
     fea.phys.ns.bdr.sel(i_outflow) = 4;
     fea.phys.ns.bdr.coef{2,end}{1,i_inflow} = s_inflow;

     fprintf( 1, '\n%s - Level %i : %g\n\n', solvers{isolv}, ilev );

% Parse and solve problem.
     if( any(isolv==[3,4])  )
       [fea.phys.ns.sfun{1:2}] = deal('sflag1');
     end
     fea = parsephys(fea);
     fea = parseprob(fea)
     fid = [];
     switch( isolv )

       case 1   % FEATool (analytic Newton Jacobian)
         jac.form  = { [1;1] [1;1] [] ;
                       [1;1] [1;1] [] ;
                       []    []    [] };
         jac.coef  = { [num2str(rho),'*ux'] [num2str(rho),'*uy'] [] ;
                       [num2str(rho),'*vx'] [num2str(rho),'*vy'] [] ;
                       []                   []                   [] };

         fid = fopen( fullfile(ldir,lfile), 'w+' );
         fea.sol.u = solvestat( fea, 'fid', fid, 'nsolve', 2, 'jac', jac );

       case 2   % Use external FEniCS solver.
         [~,file] = fileparts(lfile);
         fea = fenics( fea, 'fdir', ldir, 'fname', file, 'clean', false );

       case 3   % Use external OpenFOAM CFD solver.
         fea.sol.u = openfoam( fea, 'casedir', ldir, 'logfname', lfile, 'clean', false );

       case 4   % Use external SU2 CFD solver.
         fea.sol.u = su2( fea, 'workdir', ldir, 'logfname', lfile, 'clean', false );
     end
     [t_solv,it] = l_parse_logfile( isolv, ldir, lfile, fid );


% Calculate benchmark quantities and error.
     [cd,cl,dp,err] = l_dragliftpres( fea, rho, miu, umean, diam, i_cyl );

     nel  = size(fea.grid.c,2);
     nvt  = size(fea.grid.p,2);
     ndof = sum(fea.eqn.ndof);
     if( isolv==3 )
       ndof = 3*nel;   % OpenFOAM uses cell centered dofs.
     end
     data_i = [ data_i; [ ilev, nel, nvt, ndof, t_solv, it, cd, cl, dp, err ] ];

     delete(fullfile(ldir,lfile))
   end

   data = [ data; { data_i } ];
   if( isave )
     save 'ns3_benchmark_data' data
   end
 end


% Data processing.
 solvers  = solvers(nsolv);
 marker   = {'o','s','v','^'};
 for i=1:length(data)

   data_i = data{i};

   fprintf('\n\n%s\n|------|---------|---------|------------|---------|-----|----------------|----------------|-----------------|-----------------|----------------|\n', solvers{i} )
   fprintf('| ilev |   nel   |   nvt   |    ndof    |  t_sol  |  it |      cd_l      |      cd_v      |       cl_l      |       cl_v      |       dp       |\n' )
   fprintf('|------|---------|---------|------------|---------|-----|----------------|----------------|-----------------|-----------------|----------------|\n' )
   fprintf('| %4i | %7i | %7i | %10i | %7.1f | %3i | %14.11f | %14.11f | %15.12f | %15.12f | %14.11f |\n', data_i(:,1:end-5).' )
   fprintf('|------|---------|---------|------------|---------|-----|----------------|----------------|-----------------|-----------------|----------------|\n' )
   fprintf('| Ref. |         |         |            |         |     |  5.57953523384 |  5.57953523384 |  0.010618937712 |  0.010618937712 |  0.11752016697 |\n' )
   fprintf('|------|---------|---------|------------|---------|-----|----------------|----------------|-----------------|-----------------|----------------|\n\n\n' )

   if( iplot )
     figure(1)
     loglog( data_i(:,5), data_i(:,end-4), ['-',marker{i}], 'linewidth', 2 )
     hold on
     if( i==length(data) )
       legend( solvers{:}, 'Location', 'SouthWest' )
       grid on
       xlabel( 'CPU time [s]' )
       ylabel( 'error(c_D)' )
       title( 'Cost vs. accuracy (drag coefficient)' )
     end

     figure(2)
     loglog( data_i(:,5), data_i(:,end-2), ['-',marker{i}], 'linewidth', 2 )
     hold on
     if( i==length(data) )
       legend( solvers{:} )
       grid on
       xlabel( 'CPU time [s]' )
       ylabel( 'error(c_L)' )
       title( 'Cost vs. accuracy (lift coefficient)' )
     end

     figure(3)
     loglog( data_i(:,5), data_i(:,end), ['-',marker{i}], 'linewidth', 2 )
     hold on
     if( i==length(data) )
       legend( solvers{:} )
       grid on
       xlabel( 'CPU time [s]' )
       ylabel( 'error(\Delta p)' )
       title( 'Cost vs. accuracy (pressure difference)' )
     end
   end
 end
 if( iplot && iimg )
   figure(1)
   print -r300 -dpng ns3_benchmark_drag
   figure(2)
   print -r300 -dpng ns3_benchmark_lift
   figure(3)
   print -r300 -dpng ns3_benchmark_pres
 end


% Flow field visualization.
 if( iplot>1 )
   figure
   subplot(2,1,1)
   postplot( fea, 'surfexpr', 'sqrt(u^2+v^2)', 'arrowexpr', {'u','v'} )
   title( 'Velocity field' )
   subplot(2,1,2)
   postplot( fea, 'surfexpr', 'p' )
   title( 'Pressure' )

   if( iimg )
     print -r300 -dpng ns3_benchmark_vis
   end
 end

 if( ~nargout )
   clear data fea
 end



%