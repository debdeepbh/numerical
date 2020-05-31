%* -----------------------------------------------
%* mathspic (Perl version 1.13 Apr 26, 2010)
%* A filter program for use with PiCTeX
%* Copyright (c) 2005-2010 A Syropoulos & RWD Nickalls 
%* Command line: /usr/bin/mathspic full_proj_matlab.m
%* Input filename : full_proj_matlab.m
%* Output filename: full_proj_matlab.mt
%* Date & time: 2020/05/30    14:26:56
%* -----------------------------------------------
% generate the mesh
%script
%% *** Line 5: 
%% ***         clear all
%% ... Error: command not recognized
%% *** Line 6: 
%% ***         close all
%% ... Error: command not recognized
%delta = 0.002;	% for glass slab
%delta = 0.012;	% peridynamic horizon
%delta = 0.25;	% for the unit circle 
%% *** Line 12: 
%% ***          load('delta.mat')
%% ... Error: command not recognized
%meshsize = delta/3;
%% *** Line 15: 
%% ***          load('meshsize.mat')
%% ... Error: command not recognized
%% *** Line 18: 
%% ***          disp ('Loading matrices: Pos.mat, Vol.mat, T.mat\n')
%% ... Error: command not recognized
%% *** Line 19: 
%% ***          load('Pos.mat')
%% ... Error: command not recognized
%% *** Line 20: 
%% ***          load('Vol.mat')
%% ... Error: command not recognized
%% *** Line 21: 
%% ***          load('T.mat')
%% ... Error: command not recognized
%% *** Line 23: 
%% ***          scatter(Pos(:,1), Pos(:,2), 5, 'filled')
%% ... Error: command not recognized
%% *** Line 25: 
%% ***          disp 'Press key to continue'
%% ... Error: command not recognized
%% *** Line 26: 
%% ***          pause
%% ... Error: command not recognized
% generate the neighbor-list
%% *** Line 29: 
%% ***          [NbdArr] = gen_nbdlist(Pos, delta);
%% ... Error: command not recognized
%% *** Line 31: 
%% ***          disp('Max neighbors') 
%% ... Error: command not recognized
%% *** Line 32: 
%% ***          size(NbdArr)
%% ... Error: command not recognized
%% *** Line 33: 
%% ***          disp 'Press key to continue'
%% ... Error: command not recognized
%% *** Line 34: 
%% ***          pause
%% ... Error: command not recognized
% get the position and relative distances of the neighbors
%% *** Line 37: 
%% ***          [Nbd_xi_1, Nbd_xi_2, Nbd_xi_norm, NbdVol] = precomputation(NbdArr, Pos, Vol);
%% ... Error: command not recognized
%% *** Line 38: 
%% ***          disp 'Precomputation done.'
%% ... Error: command not recognized
% get the constant external force for the nodes
%% *** Line 41: 
%% ***          [extforce] = external_force(Pos, NbdArr);
%% ... Error: command not recognized
%% *** Line 42: 
%% ***          disp 'Press key to continue'
%% ... Error: command not recognized
%% *** Line 43: 
%% ***          pause
%% ... Error: command not recognized
% need initial condition
%% *** Line 47: 
%% ***          [NbdArr_out, u0] = simulate2(Pos, NbdArr, Vol, extforce, delta, Nbd_xi_1, Nbd_xi_2, Nbd_xi_norm);
%% ... Error: command not recognized
