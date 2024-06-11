using Plots
using DataFrames
using LaTeXStrings
using LinearAlgebra
# using GR

n::Int64 = 1000
x = range(-1,1,n)
y = range(-1,1,n)

z = zeros(Float64,n,n)
zs = zeros(Float64,n,n)

r0 = [-0.5 0 0.5]

for i=1:n
    for j=1:n
        z[i,j] = x[i]^2+y[j]^2
        zs[i,j] = 2*(x[i]-r0[1])+2*(y[j]-r0[2])-r0[3]
    end
end

surface(x, y, zs)
surface!(x,y,z)
quiver!([0],[0],[0],quiver=([-0.5],[0],[0.5]))
quiver!([-0.5],[0],[0.5],quiver=([1],[0],[-0.5]))

# plot3d(1, xlim = (-2, 2), ylim = (-2, 2),zlim = (0, 2), legend = false)
# plot3d!(x,y,z, st=[:surface , :contourf])
# l = @layout [a{0.7w} b; c{0.2h}]

Ax,Ay,Az = 2.0,3.0,7.0
Bx,By,Bz = -2.0,4.0,-1.0
A = [Ax;Ay;Az]
B = [Bx;By;Bz]
C = cross(A,B)
# Cmag = ( (Ay*Bz-By*Az)^2 + (Ax*Bz-Bx*Az)^2 + (Ax*By-Bx*Ay)^2 )^(1/2)
Cmag = norm(C)
Cnormal = 1/Cmag * C
# Cnmag = ( Cnormal[1]^2 + Cnormal[2]^2 + Cnormal[3]^2 )^(1/2)
Cnmag = norm(Cnormal)

quiver([0],[0],[0],quiver=([Ax],[Ay],[Az]))
quiver!([0],[0],[0],quiver=([Bx],[By],[Bz]))
quiver!([0],[0],[0],quiver=([C[1]],[C[2]],[C[3]]))
quiver!([0],[0],[0],quiver=([Cnormal[1]],[Cnormal[2]],[Cnormal[3]]))