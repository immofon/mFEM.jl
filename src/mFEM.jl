module mFEM

# N:int  
# Nm:int 
# P:mat 
# T:mat 
# Pb:mat, Tb:mat are 
# Nb:int 


struct Mesh
  P::AbstractArray # be an information matrix consisting of the coordinates of all mesh nodes
  Nm::Int # denotes the number of mesh nodes. Nm = size(P,1)

  T::AbstractArray # be an information matrix consisting of the global node indices of 
  N::Int # denotes the number of mesh elements. N = size(T,1)

  Integral::Function # (u::Function,n::Int)::Float64
                      # calc integral on the nth element of mesh
                      # u: (x::Tuple)::Float64
end
export Mesh

struct Basis
  Nb::Int # denote the total number of the finite element basis functions (= the number of unknowns)
  Pb::AbstractArray # the nodes corresponding to the finite element basis functions.
  Tb::AbstractArray # the nodes corresponding to the finite element basis functions.
  Nlb::Int # the number of local trial basis functions
  Dϕ::Function  # (k::Int, x::Tuple, n::Int, α::Int ) -> Real
                # return the kth deriavtive of (n,α) in mesh basis function at x
end


module Dim1
import mFEM
using ProgressMeter

# calc gauss_quad on interval [-1,1]
# u: (x::Float64)::Float64
function _gauss(u::Function)
  # c = sqrt(3/5)
  # xw = [-c 5/9; 0 8/9; c 5/9]
  
  a = 0.9061798
  b = 0.5384693
  c = 0
  wa = 0.2369269
  wb = 0.4786287
  wc = 0.5688889

  xw = [a wa; -a wa; b wb; -b wb; c wc]
  xwn,_ = size(xw)

  D = map(1:xwn) do i
    x,w = xw[i,:]
    return w * u(x)
  end
  return sum(D)
end

# calc gauss_quad on the nth element of mesh
# u: (x::Tuple)::Float64
# n index of basis element


function mesh(a,b,N)::mFEM.Mesh
  @assert N > 1
  @assert a < b

  Nm = N + 1
  h = (b - a)/N

  P = [a + (i - 1)*h for i=1:Nm]
  T = zeros(Int,N,2)
  for i=1:N
    T[i,:] = [i,i+1]
  end

  function gauss_quad(u::Function,n::Int)
    a = P[T[n,1],1]
    b = P[T[n,2],1]

    @assert a < b
    return (b-a) * _gauss(x-> u(( (b-a)*x/2 + (a+b)/2, )) ) / 2
  end

  return mFEM.Mesh(P,size(P,1),T,size(T,1),gauss_quad)
end

function linear_elem(mesh::mFEM.Mesh)::mFEM.Basis
  Pb = copy(mesh.P)
  Tb = copy(mesh.T)
  Nlb = 2

  # kth deriavtive of basis function at x
  function Dϕ(k::Int, xt::Tuple, n::Int, α::Int)
    @assert α ∈ [1,2]
    @assert k ∈ [0,1]
    @assert length(xt) == 1

    x = xt[1]

    xn1 = Pb[Tb[n,1],1]
    xn2 = Pb[Tb[n,2],1]

    if α == 1
      if k == 0
        return (xn2 - x) / (xn2 - xn1);
      elseif k == 1
        return  - 1 / (xn2 - xn1);
      else
        return 0
      end
    elseif α == 2
      if k == 0
        return (x - xn1) / (xn2 - xn1);
      elseif k == 1
        return 1 / (xn2 - xn1);
      else
        return 0
      end
    end
  end

  return mFEM.Basis(size(Pb,1),Pb,Tb,Nlb,Dϕ)
end

# - D(cDu) = f in Ω, f=g on ∂Ω
# c(x::Tuple)
function stiffness_mat(mesh::mFEM.Mesh,trial::mFEM.Basis,test::mFEM.Basis,c::Function)

  A = zeros(test.Nb,trial.Nb) # TODO sparse matrix

  @showprogress "stiffness_mat " for n = 1:mesh.N
    for α = 1:trial.Nlb, β = 1:test.Nlb

      r = mesh.Integral(n) do x::Tuple
        return c(x) * trial.Dϕ(1,x,n,α) * test.Dϕ(1,x,n,β)
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

  b = zeros(test.Nb)
  @showprogress "load_vec " for n = 1:mesh.N
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
    i,type = boundarynodes[nb,:]
    A[i,:] = zeros(nA)
    A[i,i] = 1
    b[i] = g(tuple(trial.Pb[i,:]...),type)
  end
end

function Possion(mesh::mFEM.Mesh,c::Function,f::Function,g::Function)
  trial = linear_elem(mesh)
  test = linear_elem(mesh)

  # generate boundary nodes
  boundarynodes = [1 1;trial.Nb 1]

  A = stiffness_mat(mesh,trial,test,c)
  b = load_vec(mesh,test,f)

  set_boundary_cond!(trial,boundarynodes,A,b,g)
  x = A \ b
  return x
end

function test()
  u(x) = x * cos(x)
  c(x) = exp(x)
  f(x) = - exp(x) *(cos(x) - 2sin(x) - x*cos(x) - x*sin(x))
  g(x,i) = u(x[1])

  function tuplfy1d(f::Function)::Function
    return (x::Tuple) -> f(x[1])
  end

  m = mesh(0,1,128)

  return map(u,m.P) - Possion(m,tuplfy1d(c),tuplfy1d(f),g)
end

end # Dim1

end # mFEM
