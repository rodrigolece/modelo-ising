## Implimentación del modelo de Ising con el algoritmo de Metropolis ##

module Ising

export montecarlo_mag_run, montecarlo_config_run, montecarlo_energy_run

function config_maker(L::Int)
    S=zeros(L,L)

    for i in 1:L, j in 1:L
        if rand()<0.5
            S[i,j]=-1
        else
            S[i,j]=1
        end
    end

    S
end

function spin_choose(L)
    i=rand(1:L)
    j=rand(1:L)

    i, j
end

function spin_flip!(S::Array,i,j,N=length(S),L=int(sqrt(N)))
    S[i,j]*=-1
    S
end

function periodic_energy(S,N=length(S),L=int(sqrt(N)),J=1)
    E=0

    for i in 1:L, j in 1:L
        E+=S[i,j]*(S[mod1(i-1,L),j]+S[mod1(i+1,L),j]+S[i,mod1(j-1,L)]+S[i,mod1(j+1,L)])
    end

    J*int(-E/2)

end

function ΔE(S,i,j,N=length(S),L=int(sqrt(N)),J=1)
    J*2*S[i,j]*(S[mod1(i+1,L),j]+S[mod1(i-1,L),j]+S[i,mod1(j-1,L)]+S[i,mod1(j+1,L)])
end

function acceptance(S,T,i,j,N=length(S_old),L=int(sqrt(N)),J=1)
    Δe=ΔE(S,i,j,N,L,J)
    α=exp(-(1/T)*(Δe))

    if rand()<α
        return true,Δe
    else
        return false,0
    end
end

function one_step_flip(S::Array,T,N=length(S),L=int(sqrt(N)),J=1)
    i,j=spin_choose(L)
    α,Δe=acceptance(S,T,i,j,N,L,J)
    if α==true
        spin_flip!(S,i,j,N,L)
        return Δe
    else
        return 0
    end
end

function new_energy(E_old,Δe)
    E_old+Δe
end

function montecarlo_config_run(L::Int,steps::Int,T)
    N=L*L
    S=config_maker(L)

    for i in 1:steps
        one_step_flip(S,T,N,L,1)
    end

    S
end

function montecarlo_energy_run(L::Int,steps::Int,T)
    N=L*L
    S=config_maker(L)
    E=[periodic_energy(S)]
    sizehint(E,steps)

    for i in 1:steps
        Δe=one_step_flip(S,T,N,L,1)
        push!(E,new_energy(E[i],Δe))
    end

    E
end

function montecarlo_mag_run(L::Int,steps::Int,T)
    N=L*L
    S=config_maker(L)
    M=[magnetization(S)]
    sizehint(M,steps)

    for i in 1:steps
        one_step_flip(S,T,N,L,1)
        push!(M,magnetization(S))
    end

    M
end

end