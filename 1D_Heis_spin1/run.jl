using ITensors

const N = 8          # number of sites in chain
const cutoff = 1e-10  # maximum error allowed in truncation
const maxdim = 100    # maximum bond dimension allowed

const sites = siteinds("S=1", N)   # create N spin-1 sites


# Hamiltonian of the spin-1 Heisenberg model:
function Ham(sites)
    ampo = AutoMPO()
    for j = 1:N-1
        ampo += 0.5, "S+", j, "S-", j+1
        ampo += 0.5, "S-", j, "S+", j+1
        ampo +=  1.0, "Sz", j, "Sz", j+1
    end
    return MPO(ampo, sites)
end

H = Ham(sites)

# Define a MPS
psi = randomMPS(sites, 3)


# Run DMRG algorithm
sweeps = Sweeps(5)             # define the sweeps for DMRG
setmaxdim!(sweeps, maxdim)     # set the maximum bond dimension
setcutoff!(sweeps, cutoff)     # set the maximum truncation error

energy, psi = dmrg(H, psi, sweeps)

@show energy
