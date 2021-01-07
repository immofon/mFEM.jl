module Dim2
import mFEM
using Plots

# calc gauss_quad on standrad triangle (0,0) (1,0) (0,1)
# u: (x::Float64,y::Float64)::Float64
function _gauss(u::Function)
  local points_type = [0 5 2]
  local xyz = [ 9.352701037774565E-01 3.236494811127173E-02 3.23649481112718157E-02; 
                7.612981754348137E-01 1.193509122825931E-01 0.11935091228259319E+00; 
              - 6.922209654151433E-02 5.346110482707572E-01 0.53461104827075701E+00; 
                5.933801991374367E-01 2.033099004312816E-01 0.20330990043128172E+00; 
                2.020613940682885E-01 3.989693029658558E-01 0.39896930296585570E+00; 
                5.017813831049474E-02 5.932012134282132E-01 0.35662064826129203E+00; 
                2.102201653616613E-02 8.074890031597923E-01 0.17148898030404158E+00]
  local w = [ 4.097919300803106E-02; 
              1.085536215102866E-01; 
              2.781018986881812E-03; 
              1.779689321422668E-01; 
              2.314486047444677E-01; 
              3.140226717732234E-01;
              1.242459578348437e-01]./[3,3,3,3,3,6,6]

  local ptorder3 = [1 1 2 2 3 3;
                    2 3 1 3 1 2]
  local ptorder2 = [1 2 2;
                    2 1 3]

  local v = 0.0
  for i=1:5, j=1:size(ptorder2,2)
    x = xyz[i, ptorder2[1,j]]
    y = xyz[i, ptorder2[2,j]]
    v += w[i] * u(x,y) 
  end
  for i=6:7, j=1:size(ptorder3,2)
    x = xyz[i, ptorder3[1,j]]
    y = xyz[i, ptorder3[2,j]]
    v += w[i] * u(x,y) 
  end

  return v/2 # v * measure of triangle{(0,0),(1,0),(0,1)}
end

function uniform_trianglar_mesh(x0::Vector,x1::Vector,
    N::Vector{Int})::mFEM.Mesh
  @assert size(x0,1) == size(x1,1) == 2
  @assert x0[1] < x1[1]
  @assert x0[2] < x1[2]
  @assert size(N,1) == 2
  @assert N[1] > 1 && N[2] > 1

  P = zeros((N[1]+1) * (N[2]+1),2)

  local h = (x1 - x0) ./ N
  for i in 1:(N[2]+1), j in 1:(N[1]+1)
    x,y = x0 + h .* [j-1,i-1]
    k = (i-1) * (N[1] + 1) + j
    P[k,:] = [x,y]
  end

  T = zeros(Int,2*N[1]*N[2],3)
  for i=1:N[2], j=1:N[1]
    k = (i-1) * (N[1] + 1) + j
    tk = k - (i - 1)
    T[tk,:] = [k,k+1,k+N[1]+1]
    T[tk+N[1]*N[2],:] = [k + 1 + N[1] + 1, k + N[1] + 1, k + 1]
  end


  # calc gauss_quad on the nth element of mesh
  # u: (x::Tuple)::Float64
  # n index of basis element
  function gauss_quad(u::Function,n::Int)
    local x1 = P[T[n,1],:]
    local x2 = P[T[n,2],:]
    local x3 = P[T[n,3],:]
    local O = [x2-x1 x3-x1]

    local a = x2-x1
    local b = x3-x1
    measure = sqrt(sum(a.*a)*sum(b.*b) - sum(a.*b)^2)

    return measure * _gauss((x,y)-> begin
        return u(O * [x,y] + x1)
    end)
  end
  

  return mFEM.Mesh(P,size(P,1),T,size(T,1),gauss_quad)
end

function test(n)
  tm = mFEM.Dim2.uniform_trianglar_mesh([0,0],[1,1],[n,n])
  u(x) = x[1]*x[2] *exp(x[1]^2 + x[2]^2)
  v = 0.0
  for n=1:size(tm.T,1)
    v += tm.Integral(u,n)
  end
  return abs(v - (0.5*(exp(1) - 1))^2)
end

end # Dim2

