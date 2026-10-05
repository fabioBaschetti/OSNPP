function sim = simulate_hold_closed_mc(cfg,Y)
%SIMULATE_HOLD_CLOSED_MC Keep nuclear output at the common initial level.
    Npaths=size(Y,1);
    Nt=cfg.Nt;
    P=cfg.x1_init*ones(Npaths,Nt+1);
    sim=struct('P',P,'NumSwitches',zeros(Npaths,1),'SwitchingCost',zeros(Npaths,1));
end
