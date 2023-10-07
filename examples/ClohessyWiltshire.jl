using AlgebraicControl.CMPC
using AlgebraicControl.ParaConvCat
using SCS
using Convex
using Plots

mu = 3.986e14;#standard gravtational constant[m^3/s^2]
a = big(6793137);#semi-major axis[m]
n = sqrt(mu/(a^3));#orbital rate[rad/s]
A = [0 0 0   1 0 0;
     0 0 0   0 1 0;
     0 0 0   0 0 1;
     3*(n^2) 0 0 0 2*n 0;
     0 0 0 -2*n 0 0;
     0 0 -(n^2) 0 0 0]
B = [0 0 0;
     0 0 0;
     0 0 0;
     1 0 0;
     0 1 0;
     0 0 1]
Q = 5*[1 0 0 0 0 0;
       0 1 0 0 0 0;
       0 0 1 0 0 0;
       0 0 0 1 0 0;
       0 0 0 0 1 0;
       0 0 0 0 0 1]
R = [1 0 0;
     0 1 0;
     0 0 1];#should be a 3x3 identity matrix
x₀ = [10, 10, 10, 1, 1, 1]
N = 6500
dim_x = size(A,1)
dim_u = size(R,1)

cost(uₖ,xₖ) = quadform(xₖ,Q) + quadform(uₖ,R)#change uk to quadform(uk,R)
dynamics(uₖ,xₖ) = A*xₖ + B*uₖ
#constraints(uₖ,xₖ) = [
#    uₖ <= 1, uₖ >= -1,
#    xₖ[1] <= 3, xₖ[1] >= -3,
#    xₖ[2] <= 2, xₖ[2] >= -2
#]
constraints(uₖ,xₖ)=Constraint[];#use Constraint[] to tell julia will be a vector of lists
one_step = one_step_bifunction(dim_x,dim_u,cost, constraints, dynamics)
MPC_bifunc = MPC_bifunction(one_step, N) 

us = [Variable(dim_u) for i in 1:N-1]
x_N = Variable(dim_x)

MPC_prob = to_cvx(MPC_bifunc, us, x₀, x_N);

solve!(MPC_prob, SCS.Optimizer);

function simulate(A, B, x₀, us, N)
    x = x₀
    us = evaluate.(us)
    for i in 1:N
        x = A*x + B*us[i]
    end
    return x
end

function trajectory(A, B, x₀, us, N)
     res = [x₀]
     us = evaluate.(us)
     for i in 2:N
          c = A*res[i-1] + B*us[i-1]
         push!(res, [trunc(c[i]) for i in 1:6])
     end
     return res
 end
 
 res = trajectory(A, B, x₀, us, N);

#x1=dynamics(xₖ[1],uₖ[1])#... repeat for each x
#results=[]
#push!(results,x1)
# println([us])
# print("---------------x_N---------------")
# println(x_N)
#results = [dynamics(xₖ[i],uₖ[i]) for i in 1:N-1]


#calculate using the linear time varying matrix form

#wiki eqns.
#=problem is looking for t since created and expression that uses t, but t is a changing variable not defined beforehand.
     Fixed by simply setting expression as an anonymous functon; func(....) or func = t -> ..... =# 
#establish each section of the matrix
 #xdot=🏥*x₀
xdots = []
function cl_eqs(N,n,x₀)
     for i in 1:N
          t = i
          ϕrr = [4 - 3*cos(n*t) 0 0;
          6*(sin(n*t) - n*t) 1 0;
          0 0 cos(n*t)];

          ϕrv =  [(1/n)*sin(n*t) (2/n)*(1 - cos(n*t)) 0;
          (2/n)*(cos(n*t) - 1) (1/n)*(4*sin(n*t) - 3*n*t) 0;
          0 0 (1/n)*sin(n*t)];

          ϕvr = [3*n*sin(n*t) 0 0;
          6*n*(cos(n*t) - 1) 0 0;
          0 0 -n*sin(n*t)];

          ϕvv = [cos(n*t) 2*sin(n*t) 0;
          -2*sin(n*t) (4*cos(n*t) - 3) 0;
          0 0 cos(n*t)];
          #ϕ1= setindex!(ϕ,[ϕrr ϕrv;ϕvr ϕvv])
          #xdot = [ϕrr ϕrv;ϕvr ϕvv]*x₀;
          push!(xdots,[ϕrr ϕrv;ϕvr ϕvv]*x₀)
     end
     return xdots
end
xdots = cl_eqs(N,n,x₀)
#println(xdots)
println(size(xdots,1))

plot3d(first.(xdots),getindex.(xdots,2),getindex.(xdots,3),xlabel="x",ylabel="y",zlabel="z",camera=(45,45,0))