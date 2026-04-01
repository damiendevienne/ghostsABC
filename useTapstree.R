library(JuliaCall)

julia_setup()
julia_command("using Statistics")

for (i in 1:5) {
  res <- julia_eval("mean([1,2,3,4]) + $i")
  print(res)
}