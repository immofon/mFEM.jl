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


function uniform_mesh(a,b,N)::mFEM.Mesh
  @assert N > 1
  @assert a < b

  Nm = N + 1
  h = (b - a)/N

  P = [a + (i - 1)*h for i=1:Nm]
  T = zeros(Int,N,2)
  for i=1:N
    T[i,:] = [i,i+1]
  end

  # calc gauss_quad on the nth element of mesh
  # u: (x::Tuple)::Float64
  # n index of basis element
  function gauss_quad(u::Function,n::Int)
    a = P[T[n,1],1]
    b = P[T[n,2],1]

    @assert a < b
    return (b-a) * _gauss(x-> u(( (b-a)*x/2 + (a+b)/2, )) ) / 2
  end

  return mFEM.Mesh(P,size(P,1),T,size(T,1),gauss_quad)
end


function Dϕ_reference_linear(k::Int, xt::Tuple, α::Int)
  @assert α ∈ [1,2]
  @assert k ∈ [0,1]

  local x = xt[1]

  if α == 1
    if k == 0
      return 1 - x
    elseif k == 1
      return  -1
    end
  elseif  α == 2
    if k == 0
      return x
    elseif k == 1
      return 1
    end
  end
end

function Dϕ_reference_quadratic(k::Int, xt::Tuple, α::Int)
  @assert α ∈ [1,2,3]
  @assert k ∈ [0,1,2]

  local x = xt[1]

  if α == 1
    if k == 0
      return 2x^2 - 3x + 1
    elseif k == 1
      return 4x - 3
    elseif k ==2
      return 4
    end
  elseif  α == 2
    if k == 0
      return 2x^2 - x
    elseif k == 1
      return 4x - 1
    elseif k ==2
      return 4
    end
  elseif  α == 3
    if k == 0
      return -4x^2 + 4x
    elseif k == 1
      return -8x + 4
    elseif k ==2
      return -8
    end
  end
end

function Dϕ_reference_cubic(k::Int, xt::Tuple, α::Int)
  @assert α ∈ 1:4
  @assert k ∈ 0:2

  local x = xt[1]

  if α == 1
    if k == 0
      return -4.5x^3 + 9x^2 -5.5x + 1
    elseif k == 1
      return -3*4.5x^2 + 2*9x -5.5
    elseif k == 2
      return -2*3*4.5x + 2*9
    end
  elseif α == 2
    if k == 0
      return 4.5x^3 - 4.5x^2 + x
    elseif k == 1
      return 3*4.5x^2 - 2*4.5x + 1
    elseif k == 2
      return 2*3*4.5x - 2*4.5
    end
  elseif α == 3
    if k == 0
      return 13.5x^3 - 22.5x^2 + 9x
    elseif k == 1
      return 3*13.5x^2 - 2*22.5x + 9
    elseif k == 2
      return 2*3*13.5x - 2*22.5
    end
  elseif α == 4
    if k == 0
      return -13.5x^3 + 18x^2 - 4.5x
    elseif k == 1
      return -3*13.5x^2 + 2*18x - 4.5
    elseif k == 2
      return -2*3*13.5x + 2*18
    end
  end
end

function reference2local(Dϕ_reference::Function;Pb::Array,Tb::Array,Nlb::Int)::Function
  # kth deriavtive of basis function at x
  function Dϕ(k::Int, xt::Tuple, n::Int, α::Int)
    @assert α ∈ 1:Nlb
    @assert k ∈ [0,1,2]
    @assert length(xt) == 1

    x = xt[1]

    xn1 = Pb[Tb[n,1],1]
    xn2 = Pb[Tb[n,2],1]

    x_hat = (x - xn1) / (xn2 - xn1)
    xt_hat = (x_hat,)

    if k == 0
      return Dϕ_reference(k,xt_hat,α)
    elseif k == 1
      return Dϕ_reference(k,xt_hat,α) / (xn2 - xn1)
    elseif k == 2
      @error "TODO"
    end
  end

  return Dϕ
end

function linear_elem(mesh::mFEM.Mesh)::mFEM.Basis
  Pb = copy(mesh.P)
  Tb = copy(mesh.T)
  Nlb = 2

  return mFEM.Basis(size(Pb,1),Pb,Tb,Nlb,
    reference2local(Dϕ_reference_linear;Pb=Pb,Tb=Tb,Nlb=Nlb))
end

function quadatic_elem(mesh::mFEM.Mesh)::mFEM.Basis
  Pb = zeros(Float64,size(mesh.P,1) * 2 - 1,1)
  Tb = zeros(Int,size(mesh.T,1),3)
  Nlb = 3

  for i in 1:size(mesh.P,1)
    Pb[2i - 1] = mesh.P[i]
  end

  for i in 1:(size(mesh.P,1) - 1)
    Pb[2i] = (Pb[2i-1] + Pb[2i+1]) / 2
  end
  
  for i in 1:size(mesh.T,1)
    Tb[i,1] = 2 * mesh.T[i,1] - 1
    Tb[i,2] = 2 * mesh.T[i,2] - 1
    Tb[i,3] = (Tb[i,1] + Tb[i,2]) / 2
  end

  return mFEM.Basis(size(Pb,1),Pb,Tb,Nlb,
    reference2local(Dϕ_reference_quadratic;Pb=Pb,Tb=Tb,Nlb=Nlb)
  )
end

function cubic_elem(mesh::mFEM.Mesh)::mFEM.Basis
  Pb = zeros(Float64,size(mesh.P,1) * 3 - 2,1)
  Tb = zeros(Int,size(mesh.T,1),4)
  Nlb = 4

  for i in 1:size(mesh.P,1)
    Pb[3i - 2] = mesh.P[i]
  end
  for i in 1:(size(mesh.P,1) - 1)
    Pb[3i - 1] = (2/3)* Pb[3i-2] + (1/3) * Pb[3i + 1]
    Pb[3i] = (1/3)* Pb[3i-2] + (2/3) * Pb[3i + 1]
  end

  for i in 1:size(mesh.T,1)
    Tb[i,1] = 3 * mesh.T[i,1] - 2
    Tb[i,2] = 3 * mesh.T[i,2] - 2
    Tb[i,3] = Tb[i,1] + 1
    Tb[i,4] = Tb[i,1] + 2
  end


  return mFEM.Basis(size(Pb,1),Pb,Tb,Nlb,
    reference2local(Dϕ_reference_cubic;Pb=Pb,Tb=Tb,Nlb=Nlb)
  )
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

function Possion(mesh::mFEM.Mesh,c::Function,f::Function,g::Function;
    trial::mFEM.Basis = linear_elem(mesh),test::mFEM.Basis = linear_elem(mesh))
  # generate boundary nodes
  boundarynodes = [1 1;trial.Nb 1]

  A = stiffness_mat(mesh,trial,test,c)
  b = load_vec(mesh,test,f)


  set_boundary_cond!(trial,boundarynodes,A,b,g)
  x = A \ b
  return x
end



function test(n,;finite_elem::Function=linear_elem)
  u(x) = x * cos(x)
  c(x) = exp(x)
  f(x) = - exp(x) *(cos(x) - 2sin(x) - x*cos(x) - x*sin(x))
  g(x,i) = u(x[1])

  function tuplfy1d(f::Function)::Function
    return (x::Tuple) -> f(x[1])
  end

  m = uniform_mesh(0,1,n)
  trial = finite_elem(m)
  test= finite_elem(m)

  U = Possion(m,tuplfy1d(c),tuplfy1d(f),g;trial=trial,test=test)
  return maximum(@.abs(U-map(u,trial.Pb)))
end

test_quadatic(n) = test(n,finite_elem=quadatic_elem)
test_linear(n) = test(n,finite_elem=linear_elem)
test_cubic(n) = test(n,finite_elem=cubic_elem)

end # Dim1
