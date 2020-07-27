#=
Created on Monday 01 June 2020
Last update: Wednesday 10 June 2020

@author: Michiel Stock
michielfmstock@gmail.com

Simple model based on diffusion kernel.
=#

# loading the interactions and turning them into an adjacency matrix

using DataFrames, CSV, Plots
using LinearAlgebra, StatsBase

interactions = CSV.read(joinpath("03_interaction_data", "virionette.csv"))

hosts = interactions.host_species |> unique |> sort
viruses = interactions.virus_genus |> unique |> sort

n, m = length(hosts), length(viruses)

Y = zeros(Int, n, m)

# save as df

for (h, v) in zip(interactions.host_species, interactions.virus_genus)
    hi = findfirst(x->x==h, hosts)
    vj = findfirst(x->x==v, viruses)
    Y[hi, vj] += 1
end

df = convert(DataFrame, Y)
rename!(df, Symbol.(viruses))
df.host = hosts .|> h -> replace(h, " " => "_")
CSV.write("03_interaction_data/Y.csv", df)

heatmap(Y, xlab="viruses", ylab="hosts")

"""Turns an incidence matrix `Y` in a square adjacency matrix."""
adjacencymatrix(Y::AbstractMatrix{T}) where {T<:Real} = [zeros(T, size(Y,1), size(Y,1)) Y;
                                                        Y' zeros(T, size(Y,2), size(Y,2))]

# make laplcian
laplacian(A) = Diagonal(sum(A, dims=1)[:]) - A
norm_laplacian(L) = (Diagonal(L))^(-0.5) * L * (Diagonal(L))^(-0.5)

# kernel center
centergram(K) = K .- mean(K, dims=1) .- mean(K, dims=2) .+ mean(K)

"""
    diffusion_kernel(Y; β=2.0, norm=true)

Create the diffusion kernel for an interaction matrix
"""
function diffusion_kernel(Y; β=2.0, norm=true)
    n, m = size(Y)
    # adjacency matrix
    A = adjacencymatrix(Y)
    L = laplacian(A)
    if norm
        S = diag(L) .|> (x -> x==0 ? 1.0 : x^(-0.5)) |> Diagonal
        L = S * L * S
    end
    # kernel
    K = exp(-β * L)
    return K
end

β = 2.0  # hyperparameters for random walk


K = diffusion_kernel(Y; β=2.0)

# only information for viruses
G = K[n+1:end, n+1:end]
CSV.write(joinpath("04_predictors/", "Gdiff.csv"), DataFrame(G))

Kcentered = centergram(K)

# eigenvalue decomposition for pca
# caution, eigevalues sorted from small to large
Λ, U = eigen(Symmetric(Kcentered))

plot(plot(Λ.+1e-10, yscale=:log10, title="loadings"), scatter(U[:,end-1], U[:,end-2], title="biplot"))

# embedding
Z = U * Diagonal(Λ .+ 1e-10)^0.5

Zhosts = U[1:n,:] * Diagonal(Λ .+ 1e-10)^0.5
Zviruses = U[n+1:end,:] * Diagonal(Λ .+ 1e-10)^0.5

scatter3d(Zhosts[:,end], Zhosts[:,end-1], Zhosts[:,end-2], label="hosts")
scatter3d!(Zviruses[:,end], Zviruses[:,end-1], Zviruses[:,end-2], label="viruses")
savefig("05_results/embedding3d.png")

pembh = scatter(Zhosts[:,end-1], Zhosts[:,end-2], alpha=0.2, color=:blue,
                        label=:hosts, size=(1200, 800))
annotate!([(Zhosts[i,end-1], Zhosts[i,end-2], text(hosts[i], 4)) for i in 1:n])
xlabel!("Second component")
ylabel!("Third component")
title!("Embedding hosts")

pembv = scatter(Zviruses[:,end-1], Zviruses[:,end-2], alpha=0.2, color=:orange,
                        label=:viruses, size=(1200, 800))
annotate!([(Zviruses[i,end-1], Zviruses[i,end-2], text(viruses[i], 4)) for i in 1:m])
xlabel!("Second component")
ylabel!("Third component")
title!("Embedding viruses")

plot(pembh, pembv)
savefig("05_results/embeddinghv.pdf")

# Find all likely interactions

function likelymissing(Y, K, hosts, viruses)
    n, m = size(Y)
    F = -K[1:n,end-m+1:end]
    scoredinteractions = Tuple{eltype(F),eltype(hosts),eltype(viruses)}[]
    for (j, v) in enumerate(viruses), (i, h) in enumerate(hosts)
        if Y[i,j] == 0  # negative interaction
            push!(scoredinteractions, (F[i,j], h, v))
        end
    end
    return sort!(scoredinteractions, rev=true)
end

scoredinteractions = likelymissing(Y, K, hosts, viruses)
scoredinteractions = DataFrame(scoredinteractions)
rename!(scoredinteractions, [:score, :host, :virus])

CSV.write("05_results/scoring_diffusion_kernel.csv", scoredinteractions)


# test if the diffusion kernel can predict missing interactions
# we set a fraction of the interactions to 0 and see of we can recover them

using ROC

function reconstruct(Y; β=0.1, f=0.1, norm=true)
    n, m = size(Y)
    mask = rand(n, m) .< f
    Ypert = copy(Y)
    # set all equal to 0
    Ypert[mask] .= 0
    K = diffusion_kernel(Ypert; β=β, norm=norm)
    F = K[1:n, n+1:end]
    preds = -F[mask]  # why min?
    labels = Y[mask] .> 0
    return roc(preds, labels, distances=false) |> AUC
end

trials = [reconstruct(Y, β=2.0, f=0.1, norm=true) for i in 1:100]
average_auc = mean(trials)
