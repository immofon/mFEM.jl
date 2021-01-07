module mFEM
using Printf
# N:int  
# Nm:int 
# P:mat 
# T:mat 
# Pb:mat, Tb:mat are 
# Nb:int 


struct Mesh
  P::Array # be an information matrix consisting of the coordinates of all mesh nodes
  Nm::Int # denotes the number of mesh nodes. Nm = size(P,1)

  T::Array # be an information matrix consisting of the global node indices of 
  N::Int # denotes the number of mesh elements. N = size(T,1)

  Integral::Function # (u::Function,n::Int)::Float64
                      # calc integral on the nth element of mesh
                      # u: (x::Tuple)::Float64
end
export Mesh

struct Basis
  Nb::Int # denote the total number of the finite element basis functions (= the number of unknowns)
  Pb::Array # the nodes corresponding to the finite element basis functions.
  Tb::Array # the nodes corresponding to the finite element basis functions.
  Nlb::Int # the number of local trial basis functions
  Dϕ::Function  # (k::Int, x::Tuple, n::Int, α::Int ) -> Real
                # return the kth deriavtive of (n,α) in mesh basis function at x
end


# submodules begin
include("./dim1.jl")
include("./dim2.jl")
# submodules end

# for test

# test()::Float64, error of method
function convergence_order(test::Function;Oh::Bool = true)
  N = [2^i for i=2:7]
  Error = zeros(Float64,size(N,1))
  Order = zeros(Float64,size(N,1))

  for i in 1:size(N,1)
    n = N[i]
    Error[i] = test(n)
  end

  for i in size(Error,1):-1:2
    Order[i] = Error[i-1] / Error[i]
  end

  println("   n\terror\t\torder")
  for i in 1:size(N,1)
    @printf "%4d %.10e\t%.2f" N[i] Error[i] Order[i] 
    if Oh && Order[i] > 0
      @printf "\tO(h^%.6f)" log(2,Order[i])
    end
    @printf "\n"
  end
end

end # mFEM
