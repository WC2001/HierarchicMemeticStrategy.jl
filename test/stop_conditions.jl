import hms_jl: MetaepochSummary, Deme


@testset "GlobalMetaepochLimitReached constructor validation" begin
    @test_throws ArgumentError GlobalMetaepochLimitReached(0)
    @test_throws ArgumentError GlobalMetaepochLimitReached(-1)
end


@testset "Global stopping condition meteapoch count" begin

    metaepoch_limit = 3
    global_stop_condition = GlobalMetaepochLimitReached(metaepoch_limit)

    @test !global_stop_condition([MetaepochSummary() for _ in 1:metaepoch_limit-1])
    @test  global_stop_condition([MetaepochSummary() for _ in 1:metaepoch_limit])
    @test  global_stop_condition([MetaepochSummary() for _ in 1:metaepoch_limit+1])

end


@testset "ProblemEvaluationLimitReached constructor validation" begin
    @test_throws ArgumentError ProblemEvaluationLimitReached(0)
    @test_throws ArgumentError ProblemEvaluationLimitReached(-1)
end


@testset "Global stopping condition fitness evaluation count" begin
    fitness_evaluations_limit = 1000
    global_stop_condition = ProblemEvaluationLimitReached(fitness_evaluations_limit)

    @test !global_stop_condition([MetaepochSummary(fitness_evaluation_count = fitness_evaluations_limit - 1)])
    @test global_stop_condition([MetaepochSummary(fitness_evaluation_count = fitness_evaluations_limit)])
    @test global_stop_condition([MetaepochSummary(fitness_evaluation_count = fitness_evaluations_limit + 1)])

end


@testset "MetaepochWithoutBestFitnessImprovement constructor validation" begin
    @test_throws ArgumentError MetaepochWithoutBestFitnessImprovement(0)
    @test_throws ArgumentError MetaepochWithoutBestFitnessImprovement(-1)
end

@testset "Local stopping condition metaepochs with no improvement" begin
    max_metaepochs_without_improvement = 3
    local_stop_condition = MetaepochWithoutBestFitnessImprovement(max_metaepochs_without_improvement)
    deme_root = Deme(best_fitness_values=[4.0, 3.0, 3.0, 3.0, 3.0])
    deme_stale = Deme(parent_id="test", best_fitness_values=[4.0, 3.0, 3.0, 3.0, 3.0], best_fitness=3.0)
    deme_best_fitness_live = Deme(parent_id="test", best_fitness_values=[4.0, 3.0, 3.0, 3.0], best_fitness=3.0)
    deme_best_fitness_changed = Deme(parent_id="test", best_fitness_values=[4.0, 3.0, 3.0, 3.0, 2.0], best_fitness=2.0)

    @test !local_stop_condition(deme_root, MetaepochSummary[])
    @test local_stop_condition(deme_stale, MetaepochSummary[])
    @test !local_stop_condition(deme_best_fitness_live, MetaepochSummary[])
    @test !local_stop_condition(deme_best_fitness_changed, MetaepochSummary[])

end


@testset "LocalProblemEvaluationLimitReached constructor validation" begin
    @test_throws ArgumentError LocalProblemEvaluationLimitReached(0)
    @test_throws ArgumentError LocalProblemEvaluationLimitReached(-1)
end


@testset "Local stopping condition fitness evaluation count" begin
    fitness_evaluations_limit = 1000
    local_stop_condition = LocalProblemEvaluationLimitReached(fitness_evaluations_limit)

    @test !local_stop_condition(Deme(evaluations_count=fitness_evaluations_limit-1), MetaepochSummary[])
    @test local_stop_condition(Deme(evaluations_count=fitness_evaluations_limit), MetaepochSummary[])
    @test local_stop_condition(Deme(evaluations_count=fitness_evaluations_limit+1), MetaepochSummary[])
end


@testset "Local stopping condition no active child" begin
    parent_deme = Deme(level=2)
    id = parent_deme.id
    active_child = Deme(level=3, parent_id=id)
    inactive_child = Deme(level=3, parent_id=id, is_active=false)
    local_stop_condition = AllChildrenStopped()

    metaepoch_summary_without_child = MetaepochSummary(demes=[parent_deme])
    metaepoch_summary_both_children = MetaepochSummary(demes=[parent_deme, active_child, inactive_child])
    metaepoch_summary_inactive_child = MetaepochSummary(demes=[inactive_child])

    @test !local_stop_condition(parent_deme, MetaepochSummary[])
    @test !local_stop_condition(parent_deme, [metaepoch_summary_without_child])
    @test !local_stop_condition(parent_deme, [metaepoch_summary_both_children])
    @test local_stop_condition(parent_deme, [metaepoch_summary_inactive_child])
end