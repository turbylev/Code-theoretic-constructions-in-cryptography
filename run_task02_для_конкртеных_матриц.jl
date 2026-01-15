include("isequivSSP.jl")

G1 = [
    1 0 1 1;
    0 1 0 1
]

G2 = [
    1 0 1 0;
    0 1 1 1
]

println("\nG1 (0/1):")
print01(G1)
println("\nG2 (0/1):")
print01(G2)

function check_equation(G1, G2, res)
    if res == -1
        return false
    end
    A, perm = res
    P = perm_to_P(perm)
    left  = gf2_mul(toBoolMat(A), toBoolMat(G1))
    right = gf2_mul(toBoolMat(G2), P)
    return left == right
end


run_series(sigfun, name; reps=10)

function run_series(sigfun, name; reps=10)
    println("\nЭквивалентность ($name)")

    total = 0.0
    last_res = -1

    for t in 1:reps
        elapsed = @elapsed begin
            last_res = isequivSSP(G1, G2, sigfun)
        end
        total += elapsed
    end

    avg = total / reps

    print_result(last_res)

    if last_res != -1
        ok = check_equation(G1, G2, last_res)
        println("Проверка A*G1 == G2*P: ", ok ? "ДА" : "НЕТ (ошибка)")
    end

    println("Среднее время за $reps прогонов: ", round(avg * 1000, digits=3), " мс")

    return avg, last_res
end

# Сколько прогонов
REPS = 20

t1, r1 = run_series(sig1_dim_short,      "sig1 (dim укороченного кода)"; reps=REPS)
t2, r2 = run_series(sig2_dim_hull_short, "sig2 (dim hull укороченного)"; reps=REPS)
t3, r3 = run_series(sig3_ws_hull_short,  "sig3 (спектр весов hull укороченного)"; reps=REPS)


# Итог
println("\n Итог по среднему времени ")
println("sig1: ", round(t1 * 1000, digits=3), " мс")
println("sig2: ", round(t2 * 1000, digits=3), " мс")
println("sig3: ", round(t3 * 1000, digits=3), " мс")
