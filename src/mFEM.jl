module mFEM

# N:int  
# Nm:int 
# P:mat 
# T:mat 
# Pb:mat, Tb:mat are 
# Nb:int 


struct Mesh
  N::Int # denotes the number of mesh elements
  Nm::Int # denotes the number of mesh nodes. In the one-dimentional case, Nm = N + 1
  Nb::Int # denote the total number of the finite element basis functions (= the number of unknowns)
  Nlb_trial::Int # the number of local trial basis functions
  Nlb_test::Int # the number of local test basis functions
  P::Array # be an information matrix consisting of the coordinates of all mesh nodes
  T::Array # be an information matrix consisting of the global node indices of 
  Pb::Array # the nodes corresponding to the finite element basis functions.
  Tb_trial::Array # the nodes corresponding to the finite element basis functions.
  Tb_test::Array # the nodes corresponding to the finite element basis functions.
end
export Mesh


module Dim1
import mFEM
using ProgressBars

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

  Pb = copy(P)
  Tb_trial = copy(T)
  Tb_test = copy(T)
  Nb = N + 1

  Nlb_trial = 2
  Nlb_test = 2
  
  return mFEM.Mesh(N,Nm,Nb,Nlb_trial,Nlb_test,P,T,Pb,Tb_trial,Tb_test)
end

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

# n index of basis element
function gauss_quad(u::Function,mesh::mFEM.Mesh,n::Int)
  a = mesh.P[mesh.T[n,1],1]
  b = mesh.P[mesh.T[n,2],1]

  @assert a < b
  return (b-a) * _gauss(x-> u((b-a)*x/2 + (a+b)/2)) / 2
end

# kth deriavtive of basis function at x
function Dϕ(k::Int, x::Float64, n::Int, α::Int,mesh::mFEM.Mesh)
  @assert α ∈ [1,2]
  @assert k ∈ [0,1]

  xn1 = mesh.Pb[mesh.T[n,1],1]
  xn2 = mesh.Pb[mesh.T[n,2],1]

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

# - D(cDu) = f in Ω, f=g on ∂Ω
# c(x::Float64)
function stiffness_mat(mesh::mFEM.Mesh,c::Function)

  A = zeros(mesh.Nb,mesh.Nb) # TODO sparse matrix

  for n = ProgressBar(1:mesh.N)
    for α = 1:mesh.Nlb_trial, β = 1:mesh.Nlb_test

      r = gauss_quad(mesh,n) do x::Float64
        return c(x) * Dϕ(1,x,n,α,mesh) * Dϕ(1,x,n,β,mesh)
      end

      i = mesh.Tb_test[n,β]
      j = mesh.Tb_trial[n,α]

      A[i,j] += r
    end
  end

  return A
end

# f(x::Float64)
function load_vec(mesh::mFEM.Mesh,f::Function)

  b = zeros(mesh.Nb)
  for n = 1:mesh.N
    for  β = 1:mesh.Nlb_test

      r = gauss_quad(mesh,n) do x::Float64
        return f(x) * Dϕ(0,x,n,β,mesh)
      end

      i = mesh.Tb_test[n,β]

      b[i] += r
    end
  end
  return b
end

# g(x::Float64, boundary_type::Int)
function set_boundary_cond!(mesh::mFEM.Mesh,boundarynodes::Array,
    A::Array,b::Array,g::Function)
  _,nA = size(A)
  nbn,_ = size(boundarynodes)
  for nb = 1:nbn
    i,type = boundarynodes[nb,:]
    A[i,:] = zeros(nA)
    A[i,i] = 1
    b[i] = g(mesh.Pb[i],type)
  end
end

function Possion(mesh::mFEM.Mesh,c::Function,f::Function,g::Function)
  A = stiffness_mat(mesh,c)
  b = load_vec(mesh,f)

  boundarynodes = [1 1;mesh.Nb 1]
  set_boundary_cond!(mesh,boundarynodes,A,b,g)
  x = A \ b
  return x
end

function test()
  u(x) = x * cos(x)
  c(x) = exp(x)
  f(x) = - exp(x) *(cos(x) - 2sin(x) - x*cos(x) - x*sin(x))
  g(x,i) = u(x)

  m = mesh(0,1,10000)

  return map(u,m.P) - Possion(m,c,f,g)
end

end # Dim1

end # mFEM
