using Symbolics
using LinearAlgebra

@variables u v c a R

@variables u v t

∂ᵤ = Differential(u)
∂ᵥ = Differential(v)

∂ᵤ = Differential(u)
∂ᵥ = Differential(v)
∂ₜ = Differential(t)

#r = [u, v , sin(u + v)]
r = [u, v , sin(u) + sin(v)]

#r = [(c+a*cos(v))*cos(u), (c+a*cos(v))*sin(u), a*sin(v)] #Torus

#r = [R*sin(u)*cos(v), R*sin(u)*sin(v), R*cos(u)] #Sphere

#r = [u, v, t, sin(u + v) + sin(t)]


rᵤ = [expand_derivatives(∂ᵤ(r[1])), expand_derivatives(∂ᵤ(r[2])), expand_derivatives(∂ᵤ(r[3]))]

rᵥ = [expand_derivatives(∂ᵥ(r[1])), expand_derivatives(∂ᵥ(r[2])), expand_derivatives(∂ᵥ(r[3]))]

rₜ = [expand_derivatives(∂ₜ(r[1])), expand_derivatives(∂ₜ(r[2])), expand_derivatives(∂ₜ(r[3])), expand_derivatives(∂ₜ(r[4]))]

rᵣᵤₓᵣᵥ = cross(rᵤ,rᵥ)

rₙₒᵣₘ = norm(rᵣᵤₓᵣᵥ)

"""
gᵢⱼ = [1 + (cos(u+v))^2 (cos(u+v))^2
        (cos(u+v))^2     1 + (cos(u+v))^2]
"""
gᵢⱼ = [1 + (cos(u))^2 cos(u)*cos(v)
        cos(u)*cos(v) 1 + (cos(v))^2]
"""
#2-Torus
gᵢⱼ = [(c + a*cos(v))^2 0
        0     a^2]

#2-Sphere
gᵢⱼ = [R^2 0
        0  R^2*(sin(u))^2]


gᵢⱼ = [1+(cos(u+v))^2 (cos(u+v))^2 cos(u+v)*cos(t)
       (cos(u+v))^2 1+(cos(u+v))^2 cos(u+v)*cos(t)
       cos(u+v)*cos(t) cos(u+v)*cos(t) 1+(cos(t))^2]
"""

"""
gⁱʲ = (1/(1 + 2(cos(u+v))^2))*[1 + (cos(u+v))^2 -(cos(u+v))^2
        -(cos(u+v))^2     1 + (cos(u+v))^2]
"""
gⁱʲ = (1/(1+(cos(u))^2+(cos(v))^2))*[1 + (cos(u))^2 -cos(u)*cos(v)
                                    -cos(u)*cos(v) 1 + (cos(v))^2]
"""

gⁱʲ = [1/(c+a*cos(v))^2 0
       0 1/(a^2)]


gⁱʲ = [1/R^2 0
        0    1/(R^2*(sin(u))^2)]
"""

gⁱʲ

gᵢⱼᵤ = [expand_derivatives(∂ᵤ(gᵢⱼ[1,1])) expand_derivatives(∂ᵤ(gᵢⱼ[1,2]))
         expand_derivatives(∂ᵤ(gᵢⱼ[2,1])) expand_derivatives(∂ᵤ(gᵢⱼ[2,2]))]

gᵢⱼᵥ = [expand_derivatives(∂ᵥ(gᵢⱼ[1,1])) expand_derivatives(∂ᵥ(gᵢⱼ[1,2]))
         expand_derivatives(∂ᵥ(gᵢⱼ[2,1])) expand_derivatives(∂ᵥ(gᵢⱼ[2,2]))]

N̂ = rᵣᵤₓᵣᵥ / rₙₒᵣₘ

N̂ = [-sin(u)*cos(v), -sin(u)*sin(v), -cos(u)]

N̂ᵤ = [expand_derivatives(∂ᵤ(N̂[1])), expand_derivatives(∂ᵤ(N̂[2])), expand_derivatives(∂ᵤ(N̂[3]))]

N̂ᵥ = [expand_derivatives(∂ᵥ(N̂[1])), expand_derivatives(∂ᵥ(N̂[2])), expand_derivatives(∂ᵥ(N̂[3]))]

B = [dot(-N̂ᵤ, rᵤ) dot(-N̂ᵥ, rᵤ)
     dot(-N̂ᵤ, rᵥ) dot(-N̂ᵥ, rᵥ)]

B = [R 0
     0 (sin(u))^2*R]

C = [dot(rᵤ, rᵤ) dot(rᵥ, rᵤ)
     dot(rᵤ, rᵥ) dot(rᵥ, rᵥ)]

C = [R^2 0
     0   (sin(u))^2*R^2]

C⁻¹ = gⁱʲ

shape_matrix = C⁻¹ * B

K = det(shape_matrix)

H = (1/2) * tr(shape_matrix)

gⁱʲ

gᵢⱼᵤ

gᵢⱼᵥ

Γᵘᵤᵤ = (1/2) * (gⁱʲ[1,1] * (gᵢⱼᵤ[1,1] + gᵢⱼᵤ[1,1] - gᵢⱼᵤ[1,1]) + gⁱʲ[1,2] * (gᵢⱼᵤ[2,1] + gᵢⱼᵤ[2,1] - gᵢⱼᵥ[1,1]))

Γᵘᵤᵥ = (1/2) * (gⁱʲ[1,1] * (gᵢⱼᵥ[1,1] + gᵢⱼᵤ[1,2] - gᵢⱼᵤ[1,2]) + gⁱʲ[1,2] * (gᵢⱼᵥ[2,1] + gᵢⱼᵤ[2,2] - gᵢⱼᵥ[1,2]))

Γᵘᵥᵤ = (1/2) * (gⁱʲ[1,1] * (gᵢⱼᵤ[1,2] + gᵢⱼᵥ[1,1] - gᵢⱼᵤ[2,1]) + gⁱʲ[1,2] * (gᵢⱼᵤ[2,2] + gᵢⱼᵥ[2,1] - gᵢⱼᵥ[2,1]))

Γᵘᵥᵥ = (1/2) * (gⁱʲ[1,1] * (gᵢⱼᵥ[1,2] + gᵢⱼᵥ[1,2] - gᵢⱼᵤ[2,2]) + gⁱʲ[1,2] * (gᵢⱼᵥ[2,2] + gᵢⱼᵥ[2,2] - gᵢⱼᵥ[2,2]))

Γᵛᵤᵤ = (1/2) * (gⁱʲ[2,1] * (gᵢⱼᵤ[2,1] + gᵢⱼᵤ[2,1] - gᵢⱼᵥ[1,1]) + gⁱʲ[2,2] * (gᵢⱼᵤ[2,1] + gᵢⱼᵤ[2,1] - gᵢⱼᵥ[1,1]))

Γᵛᵤᵥ = (1/2) * (gⁱʲ[2,1] * (gᵢⱼᵤ[1,2] + gᵢⱼᵥ[1,1] - gᵢⱼᵤ[1,2]) + gⁱʲ[2,2] * (gᵢⱼᵤ[2,2] + gᵢⱼᵥ[2,1] - gᵢⱼᵥ[1,2]))

Γᵛᵥᵤ = (1/2) * (gⁱʲ[2,1] * (gᵢⱼᵤ[1,2] + gᵢⱼᵥ[1,1] - gᵢⱼᵤ[2,1]) + gⁱʲ[2,2] * (gᵢⱼᵤ[2,2] + gᵢⱼᵥ[2,1] - gᵢⱼᵥ[2,1]))

Γᵛᵥᵥ = (1/2) * (gⁱʲ[2,1] * (gᵢⱼᵥ[1,2] + gᵢⱼᵥ[1,2] - gᵢⱼᵤ[2,2]) + gⁱʲ[2,2] * (gᵢⱼᵥ[2,2] + gᵢⱼᵥ[2,2] - gᵢⱼᵥ[2,2]))

Γᵘᵤᵤᵤ = expand_derivatives(∂ᵤ(Γᵘᵤᵤ))
Γᵘᵤᵤᵥ = expand_derivatives(∂ᵥ(Γᵘᵤᵤ))

Γᵘᵤᵥᵤ = expand_derivatives(∂ᵤ(Γᵘᵤᵥ))
Γᵘᵤᵥᵥ = expand_derivatives(∂ᵥ(Γᵘᵤᵥ))

Γᵘᵥᵤᵤ = expand_derivatives(∂ᵤ(Γᵘᵥᵤ))
Γᵘᵥᵤᵥ = expand_derivatives(∂ᵥ(Γᵘᵥᵤ))

Γᵘᵥᵥᵤ = expand_derivatives(∂ᵤ(Γᵘᵥᵥ))
Γᵘᵥᵥᵥ = expand_derivatives(∂ᵥ(Γᵘᵥᵥ))

Γᵛᵤᵤᵤ = expand_derivatives(∂ᵤ(Γᵛᵤᵤ))
Γᵛᵤᵤᵥ = expand_derivatives(∂ᵥ(Γᵛᵤᵤ))

Γᵛᵤᵥᵤ = expand_derivatives(∂ᵤ(Γᵛᵤᵥ))
Γᵛᵤᵥᵥ = expand_derivatives(∂ᵥ(Γᵛᵤᵥ))

Γᵛᵥᵤᵤ = expand_derivatives(∂ᵤ(Γᵛᵥᵤ))
Γᵛᵥᵤᵥ = expand_derivatives(∂ᵥ(Γᵛᵥᵤ))

Γᵛᵥᵥᵤ = expand_derivatives(∂ᵤ(Γᵛᵥᵥ))
Γᵛᵥᵥᵥ = expand_derivatives(∂ᵥ(Γᵛᵥᵥ))

Rᵘᵤᵤᵤ = Γᵘᵤᵤᵤ - Γᵘᵤᵤᵤ - Γᵘᵤᵤ*Γᵘᵤᵤ - Γᵘᵥᵤ*Γᵛᵤᵤ  + Γᵘᵤᵤ*Γᵘᵤᵤ + Γᵘᵥᵤ*Γᵛᵤᵤ

Rᵘᵤᵤᵥ = Γᵘᵤᵥᵤ - Γᵘᵤᵤᵥ - Γᵘᵤᵥ*Γᵘᵤᵤ - Γᵘᵥᵥ*Γᵛᵤᵤ + Γᵘᵤᵤ*Γᵘᵤᵥ + Γᵘᵥᵤ*Γᵛᵤᵥ
Rᵘᵤᵥᵤ = -Rᵘᵤᵤᵥ

Rᵘᵤᵥᵥ = Γᵘᵤᵥᵥ - Γᵘᵤᵥᵥ - Γᵘᵤᵥ*Γᵘᵤᵥ - Γᵘᵥᵥ*Γᵛᵤᵥ + Γᵘᵤᵥ*Γᵘᵤᵥ + Γᵘᵥᵥ*Γᵛᵤᵥ

Rᵘᵥᵤᵤ = Γᵘᵥᵤᵤ - Γᵘᵥᵤᵤ - Γᵘᵤᵤ*Γᵘᵥᵤ - Γᵘᵥᵤ*Γᵛᵥᵤ + Γᵘᵤᵤ*Γᵘᵥᵤ + Γᵘᵥᵤ*Γᵛᵥᵤ

Rᵘᵥᵤᵥ = Γᵘᵥᵥᵤ - Γᵘᵥᵤᵥ - Γᵘᵤᵥ*Γᵘᵥᵤ - Γᵘᵥᵥ*Γᵛᵥᵤ + Γᵘᵤᵤ*Γᵘᵥᵥ + Γᵘᵥᵤ*Γᵛᵥᵥ
Rᵘᵥᵥᵤ = -Rᵘᵥᵤᵥ

Rᵘᵥᵥᵥ = Γᵘᵥᵥᵥ - Γᵘᵥᵥᵥ - Γᵘᵤᵥ*Γᵘᵥᵥ - Γᵘᵥᵥ*Γᵛᵥᵥ + Γᵘᵤᵥ*Γᵘᵥᵥ + Γᵘᵥᵥ*Γᵛᵥᵥ

Rᵛᵤᵤᵤ = Γᵛᵤᵤᵤ - Γᵛᵤᵤᵤ - Γᵛᵤᵤ*Γᵘᵤᵤ - Γᵛᵥᵤ*Γᵛᵤᵤ  + Γᵛᵤᵤ*Γᵘᵤᵤ + Γᵛᵥᵤ*Γᵛᵤᵤ

Rᵛᵤᵤᵥ = Γᵛᵤᵥᵤ - Γᵛᵤᵤᵥ - Γᵛᵤᵥ*Γᵘᵤᵤ - Γᵛᵥᵥ*Γᵛᵤᵤ + Γᵛᵤᵤ*Γᵘᵤᵥ + Γᵛᵥᵤ*Γᵛᵤᵥ
Rᵛᵤᵥᵤ = -Rᵛᵤᵤᵥ

Rᵛᵤᵥᵥ = Γᵛᵤᵥᵥ - Γᵛᵤᵥᵥ - Γᵛᵤᵥ*Γᵘᵤᵥ - Γᵛᵥᵥ*Γᵛᵤᵥ + Γᵛᵤᵥ*Γᵘᵤᵥ + Γᵛᵥᵥ*Γᵛᵤᵥ

Rᵛᵥᵤᵤ = Γᵛᵥᵤᵤ - Γᵛᵥᵤᵤ - Γᵛᵤᵤ*Γᵘᵥᵤ - Γᵛᵥᵤ*Γᵛᵥᵤ + Γᵛᵤᵤ*Γᵘᵥᵤ + Γᵛᵥᵤ*Γᵛᵥᵤ

Rᵛᵥᵤᵥ = Γᵛᵥᵥᵤ - Γᵛᵥᵤᵥ - Γᵛᵤᵥ*Γᵘᵥᵤ - Γᵛᵥᵥ*Γᵛᵥᵤ + Γᵛᵤᵤ*Γᵘᵥᵥ + Γᵛᵥᵤ*Γᵛᵥᵥ
Rᵛᵥᵥᵤ = -Rᵛᵥᵤᵥ

Rᵛᵥᵥᵥ = Γᵛᵥᵥᵥ - Γᵛᵥᵥᵥ - Γᵛᵤᵥ*Γᵘᵥᵥ - Γᵛᵥᵥ*Γᵛᵥᵥ + Γᵛᵤᵥ*Γᵘᵥᵥ + Γᵛᵥᵥ*Γᵛᵥᵥ

Rᵤᵤ = Rᵘᵤᵤᵤ + Rᵛᵤᵥᵤ

Rᵤᵥ = Rᵘᵤᵤᵥ + Rᵛᵤᵥᵥ

Rᵥᵤ = Rᵘᵥᵤᵤ + Rᵛᵥᵥᵤ

Rᵥᵥ = Rᵘᵥᵤᵥ + Rᵛᵥᵥᵥ

R = gⁱʲ[1,1]*Rᵤᵤ + gⁱʲ[1,2]*Rᵤᵥ + gⁱʲ[2,1]*Rᵥᵤ + gⁱʲ[2,2]*Rᵥᵥ

gⁱʲ[1,1]

Rᵤᵤ

"""
gⁱʲ = [1/(c+a*cos(v))^2 0
       0 1/(a^2)]

gᵢⱼᵤ = [0 0
        0 0]

gᵢⱼᵥ = [-2*a*sin(v)*(c+a*cos(v)) 0
        0                        0]
"""


