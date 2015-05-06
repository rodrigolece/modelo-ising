## Implimentación del modelo de Ising con el algoritmo de Metropolis ##

module Ising

export MicroEstado, edo_aleatorio, simulacion_montecarlo, montecarlo_energia, monte
export voltea_espin!, energia_total, energia_ij, propone_cambio, paso_montecarlo
import Base.show

type MicroEstado
    σ::Array{Int,2}
	#Vamos a suponer que todas las configuraciones son cuadradas
    L::Int
end

show(io::IO, m::MicroEstado) = print(io, m.σ)

function edo_aleatorio(L::Int)
    σ = ones(Int, (L,L))
    for i in 1:L^2
        if rand() <= 0.5
			σ[i] = -1
        end
	end
    MicroEstado(σ,L)
end

function voltea_espin!(m::MicroEstado, i::Int, j::Int)
    m.σ[i,j] *= -1
end

function energia_total(m::MicroEstado)
	out = 0.
	L = m.L
    for i in 1:L, j in L
		out -= m.σ[i,j]*(m.σ[mod1(i-1,L),j] + m.σ[mod1(i+1,L),j] + m.σ[i,mod1(j-1,L)] + m.σ[i,mod1(j+1,L)])
    end
    out/2
end

function energia_ij(m::MicroEstado, i::Int, j::Int)
	L = m.L
	-0.5 * m.σ[i,j]*(m.σ[mod1(i-1,L),j] + m.σ[mod1(i+1,L),j] + m.σ[i,mod1(j-1,L)] + m.σ[i,mod1(j+1,L)])
end

function propone_cambio(m::MicroEstado, β::Float64)
    i, j = rand(1:m.L), rand(1:m.L)  #Es más rápido que rand(1:m.L, 2)
	ΔE = -2*energia_ij(m, i, j)

	ΔE, i, j
end

function paso_montecarlo(m::MicroEstado, β::Float64)
	ΔE, i, j = propone_cambio(m, β)

	#El parámetro 	de aceptación
	α = min(1., e^(-β*ΔE))

    if rand() < α
        voltea_espin!(m, i, j)
        return ΔE
    else
        return 0.
    end
end

function simulacion_montecarlo(L::Int, T::Float64, num_pasos::Int)
	β = 1/T
	m = edo_aleatorio(L)

	for i in 1:num_pasos
		paso_montecarlo(m,β)
	end

	m
end

function montecarlo_energia(L::Int, T::Float64, num_pasos::Int)
	β = 1/T
	m = edo_aleatorio(L)

	out = [energia_total(m)]
    sizehint(out, num_pasos)

    for i in 1:num_pasos-1
        ΔE = paso_montecarlo(m, β)
        push!(out, out[i] + ΔE)
    end

    out
end

magnetizacion(m::MicroEstado) = sum(m.σ)

function montecarlo_magnetizacion(L::Int, T::Float64, num_pasos::Int)
    β = 1/T
	m = edo_aleatorio(L)

    out = [magnetizacion(m)]
    sizehint(out, num_pasos)

    for i in 1:num_pasos-1
        paso_montecarlo(m, β)
        push!(out, magnetizacion(m))
    end

    out
end



# function microestados(n :: Int64, m :: Int64)
#     N=n*m
#     if N>16 return 0 end
#     N2=big(2)^N
#     out=zeros(N2,N)
#     for i in 1:N2-1
#         bini=bin(i)
#         out[i+1,N-length(bini)+1:N]=(int(split(bini,"")).*2)-1
#     end
#     out
# end

# function particion_T(T :: Float64, configuraciones :: Array{Float64,2}, n :: Int64, m :: Int64)
#     out=0.0
#     for i in length(configuraciones[:,1])
#         out+=e^(-energia_total(configuraciones[i,:],n,m)/T)
#     end
#     out
# end

end
