## Implimentación del modelo de Ising con el algoritmo de Metropolis-Hastings ##

module Ising

export MicroEstado, edo_aleatorio, energia_total, magnetizacion_total, simulacion_montecarlo
export microEstados_montecarlo
import Base.show

type MicroEstado
    σ::Array{Int,2}
	# Vamos a suponer que todas las configuraciones son cuadradas
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
    for i in 1:m.L, j in m.L
		out += energia_ij(m,i,j)
    end
    out/2
end

function energia_ij(m::MicroEstado, i::Int, j::Int)
	L = m.L
	-m.σ[i,j]*(m.σ[mod1(i-1,L),j] + m.σ[mod1(i+1,L),j] + m.σ[i,mod1(j-1,L)] + m.σ[i,mod1(j+1,L)])
end

function propone_cambio(m::MicroEstado, β::Float64)
    i, j = rand(1:m.L), rand(1:m.L)  # Es más rápido que rand(1:m.L, 2)
	ΔE = -2*energia_ij(m, i, j)

	ΔE, i, j
end

function paso_montecarlo!(m::MicroEstado, β::Float64)
	aceptado = false

	while aceptado == false
		ΔE, i, j = propone_cambio(m, β)

		# El parámetro de aceptación
		α = min(1., e^(-β*ΔE))

		if rand() < α
			aceptado = true
			ΔM = -2*m.σ[i,j]
			voltea_espin!(m, i, j)
			return ΔE, ΔM
		end
    end
end

# function simulacion_montecarlo(L::Int, T, num_pasos::Int)
# 	β = 1/T
# 	m = edo_aleatorio(L)

# 	for i in 1:num_pasos-1
# 		paso_montecarlo(m,β)
# 	end

# 	m
# end

magnetizacion_total(m::MicroEstado) = sum(m.σ)

function simulacion_montecarlo(L::Int, T, num_pasos::Int)
	β = 1/T
	m = edo_aleatorio(L)

	ener = Array(Float64, num_pasos)
	ener[1] = energia_total(m)
	mag = Array(Float64, num_pasos)
	mag[1] = magnetizacion_total(m)

	for i in 1:num_pasos-1
		ΔE, ΔM = paso_montecarlo!(m, β)
		ener[i+1] = ener[i] + ΔE
		mag[i+1] = mag[i] + ΔM
	end

	ener, mag
end

function microEstados_montecarlo(L::Int, T, num_pasos::Int)
	β = 1/T
	m = edo_aleatorio(L)

	out = Array{Int,2}[copy(m.σ)]
	sizehint(out, num_pasos)

	for i in 1:num_pasos-1
		paso_montecarlo!(m, β)
		push!(out, copy(m.σ))
	end

	out
end

end
