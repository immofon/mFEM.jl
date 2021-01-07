module mFEM
using Printf
using ProgressMeter

struct Mesh
  P::Array # be an information matrix consisting of the coordinates of all mesh nodes
  T::Array # be an information matrix consisting of the global node indices 
  BNode::Array # be an information matrix consisting to the node indices of boundary nodes
                # e.g. [0,2,3]
  BEage::Array # e.g. [0 1;3 4;2 5], with counterclockwise 
  BSurface::Array # e.g. [0 1 3; 2 3 5] with counterclockwise
  Integral::Function # (u::Function,n::Int)::Float64
                      # calc integral on the nth element of mesh
                      # u: (x::Tuple)::Float64
end
export Mesh

export Basis
struct Basis
  Nlb::Int # the number of local trial basis functions
  Pb::Array # the nodes corresponding to the finite element basis functions.
  Tb::Array # the nodes corresponding to the finite element basis functions.
  BNode::Array # be an information matrix consisting to the node indices of boundary nodes
                # e.g. [0,2,3]
  BEage::Array # e.g. [0 1;3 4;2 5], with counterclockwise 
  BSurface::Array # e.g. [0 1 3; 2 3 5] with counterclockwise
  Dϕ::Function  # (k::Int, x::Tuple, n::Int, α::Int ) -> Real
                # return the kth deriavtive of (n,α) in mesh basis function at x
end

export stiffness_mat,load_vec,set_boundary_cond!
# - D(cDu) = f in Ω, f=g on ∂Ω
# c(x::Tuple)
function stiffness_mat(mesh::mFEM.Mesh,trial::mFEM.Basis,test::mFEM.Basis,c::Function)

  A = zeros(size(test.Pb,1),size(trial.Pb,1)) # TODO sparse matrix

  @showprogress "stiffness_mat " for n = 1:size(mesh.T,1)
    for α = 1:trial.Nlb, β = 1:test.Nlb

      r = mesh.Integral(n) do x::Tuple
        return c(x) * trial.Dϕ(1,x,n,α) .* test.Dϕ(1,x,n,β)
      end

      i = test.Tb[n,β]
      j = trial.Tb[n,α]

      A[i,j] += r
    end
  end

  return A
end

# f(x::Tuple)
function load_vec(mesh::mFEM.Mesh,test::mFEM.Basis,f::Function)

  b = zeros(size(test.Pb,1))
  @showprogress "load_vec " for n = 1:size(mesh.T,1)
    for  β = 1:test.Nlb

      r = mesh.Integral(n) do x::Tuple
        return f(x) * test.Dϕ(0,x,n,β)
      end

      i = test.Tb[n,β]

      b[i] += r
    end
  end
  return b
end

# g(x::Tuple, boundary_type::Int)
function set_boundary_cond!(trial::mFEM.Basis,boundarynodes::AbstractArray,
    A::AbstractArray,b::AbstractArray,g::Function)
  nA = size(A,2)
  nbn = size(boundarynodes,1)
  for nb = 1:nbn
    type,i = boundarynodes[nb,:]
    if type == -1 # Dirichlet
      A[i,:] = zeros(nA)
      A[i,i] = 1
      b[i] = g(tuple(trial.Pb[i,:]...),type)
    else
      error("unsupported type")
    end
  end
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
