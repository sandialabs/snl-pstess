if k >= sum(k_inc(1:3))+1
       switch_cr = find(kscr_sw == k);
       if ~isempty(switch_cr)
          if switch_cr == 1
            bus_pf2 = bus_scr1;
            Y_gpf2 = Y_gscr1;
            Y_gncpf2 = Y_gncscr1;
            Y_ncgpf2 = Y_ncgscr1;
            Y_ncpf2 = Y_ncscr1;
            V_rgpf2 = V_rgscr1;
            V_rncpf2 = V_rncscr1;
            boscr2 = boscr1;
            bus_intpf2 = bus_intscr1;
          elseif switch_cr == 2
            bus_pf2 = bus_scr2;
            Y_gpf2 = Y_gscr2;
            Y_gncpf2 = Y_gncscr2;
            Y_ncgpf2 = Y_ncgscr2;
            Y_ncpf2 = Y_ncscr2;
            V_rgpf2 = V_rgscr2;
            V_rncpf2 = V_rncscr2;
            boscr2 = boscr2;
            bus_intpf2 = bus_intscr2;
          elseif switch_cr == 3
            bus_pf2 = bus_scr3;
            Y_gpf2 = Y_gscr3;
            Y_gncpf2 = Y_gncscr3;
            Y_ncgpf2 = Y_ncgscr3;
            Y_ncpf2 = Y_ncscr3;
            V_rgpf2 = V_rgscr3;
            V_rncpf2 = V_rncscr3;
            boscr2 = boscr3;
            bus_intpf2 = bus_intscr3;
          elseif switch_cr == 4
            bus_pf2 = bus_scr4;
            Y_gpf2 = Y_gscr4;
            Y_gncpf2 = Y_gncscr4;
            Y_ncgpf2 = Y_ncgscr4;
            Y_ncpf2 = Y_ncscr4;
            V_rgpf2 = V_rgscr4;
            V_rncpf2 = V_rncscr4;
            boscr2 = boscr4;
            bus_intpf2 = bus_intscr4;
          end
       end