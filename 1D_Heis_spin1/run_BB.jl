using ITensors

const N = 16          # number of sites in chain
const cutoff = 1e-10  # maximum error allowed in truncation
const maxdim = 100    # maximum bond dimension allowed

const sites = siteinds("S=1", N)   # create N spin-1 sites

beta=1
# Hamiltonian of the spin-1 Heisenberg model:
function Ham(sites)
    ampo = AutoMPO()
    for j = 1:N-1
        # Bilinear term
        ampo += 0.5, "S+", j, "S-", j+1
        ampo += 0.5, "S-", j, "S+", j+1
        ampo +=  1.0, "Sz", j, "Sz", j+1

        #Biquadratic term
        ampo += beta, "Sz", j, "Sz", j+1, "Sz", j, "Sz", j+1
        ampo += beta*0.5, "Sz", j, "Sz", j+1, "S+", j, "S-", j+1
        ampo += beta*0.5, "Sz", j, "Sz", j+1, "S-", j, "S+", j+1

        ampo += beta*0.5, "S+", j, "S-", j+1, "Sz", j, "Sz", j+1
        ampo += beta*0.25, "S+", j, "S-", j+1, "S-", j, "S+", j+1
        ampo += beta*0.25, "S+", j, "S-", j+1, "S+", j, "S-", j+1
        
        ampo += beta*0.5, "S-", j, "S+", j+1, "Sz", j, "Sz", j+1
        ampo += beta*0.25, "S-", j, "S+", j+1, "S-", j, "S+", j+1
        ampo += beta*0.25, "S-", j, "S+", j+1, "S+", j, "S-", j+1
    end
    return MPO(ampo, sites)
end

H = Ham(sites)

# if beta==1/3
#     H = H+2/3
# elseif beta==2
#     H = H-12
# elseif beta==1
#     H= H-4
# end

# Define a MPS
psi = randomMPS(sites, 3)


# Run DMRG algorithm
sweeps = Sweeps(50)             # define the sweeps for DMRG
setmaxdim!(sweeps, maxdim)     # set the maximum bond dimension
setcutoff!(sweeps, cutoff)     # set the maximum truncation error

energy, psi = dmrg(H, psi, sweeps)

if beta==1/3
    @show energy+2/3
elseif beta==2
    @show energy-12
elseif beta==1
    @show energy-4
end
#@show energy
