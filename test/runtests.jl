using Kuramoto, Test
Km = Kuramoto

using Random, Distributions, DynamicalSystems

ω_dist = Normal(0, 1)
Kc = Km.critical(ω_dist)

@test Kc ≈ 1.5957691216057306

Ks = [0, Kc / 2, Kc, 3Kc / 2, 2Kc]
rs = [
    0.0642670918200937,
    0.16226278955464907,
    0.448649488627585,
    0.8411925130146383,
    0.9541400637196635
]

Random.seed!(1)

N = 10^2
T = 10^3
dt = 0.05

for (K, r) in zip(Ks, rs)
    println("Simulating with N=$N, T=$T and K=$K...")
    test_system = Km.system(N, ω_dist, K=K)
    trj = trajectory(test_system, T, dt=dt)
    @test Km.order(trj)[end] ≈ r
end
