% What I can do now is to create a wrapper function solver for each 
% student submission and test the solver function here
classdef autograder < matlab.unittest.TestCase 
    properties
        tol = 1e-5
        C_random = (rand(1, 3)-0.5)*20
        coeff
        solver
    end
    
    methods(TestClassSetup)
        function add_solvers(testCase)
            testCase.solver = generate_solvers();
        end
    end
    methods(Test)
        function test_random_case(testCase)
            testCase.coeff = [1, testCase.C_random, -1];
            fields = fieldnames(testCase.solver);
            for i = 1:length(fields)
                fprintf('Testing solver : %s', fields{i});
                s = testCase.solver.(fields{i});
                rts = s(testCase.C_random);
                for j = 1:length(rts)
                    testCase.verifyEqual(polyval(testCase.coeff, rts(j)), ...
                        0.0, 'AbsTol', testCase.tol);
                end
            end
        end
    end
end