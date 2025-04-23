import numpy as np
import pdb
import argparse
import subprocess # For executing c++ executable
import pandas as pd
from timeit import default_timer as timer

###############################################################
################### Util Functions Below ######################

def convertPIs(aString):
    """ Input: A comma seperated string like "pi/2,pi/4,pi/2,pi/4,pi/2,"
            or 1.2,4,5.3 etc
    Output: string replacing the pis 1.57079,...,...
    """
    if aString[-1] == ",":  # Remove training comma if there is one
        aString = aString[:-1]
    aString = aString.replace("pi", "3.141592") # Replace pi with 3.14... if needed
    vecOfStrings = aString.split(",")
    ans = []
    for anExpression in vecOfStrings:
        ans.append(str(eval(anExpression))) # Evaluate expressions if needed
    return ans

###############################################################
################### Main Functions Below ######################


def graderMain(executablePath, gradingCSV):
    problems = [[
        "./map2.txt",
        "1.1066466592969404,0.9723013030643942,0.4952398983544146,2.7652695620174006,2.0395622038882006,2.6214497067815428,0.7867338202664255,0.2636248905397965",
        "0.6466152543377316,2.526447753225196,0.8089663287312084,2.75260909193042,2.520836055574241,0.38273522398036697,0.882453764993128,0.4109044557604403"
    ], [
        "./map2.txt",
        "1.8431104532359168,0.33563942536758407,2.7973887009956524,0.35783682421496366,1.8995954921019715,2.655846293781794,1.6579980294170018,0.5757484019804435",
        "0.28052659794416446,2.0043176674035776,2.6212711431331592,1.9426172464798803,0.18682440358961935,1.557193503240161,2.6597079814300226,0.804962434925884"
    ], [
        "./map2.txt",
        "1.2953818367655228,2.58018416204552,0.5123162539397882,2.4884060599827875,0.278032475269281,2.345344975104387,0.29404410522675095,0.3530815644619991",
        "0.2401068244236702,1.692168449298867,2.0509271850506323,0.5100714156510705,2.738218698011262,2.1966322241228426,1.240080805920904,0.213566944739119"
    ]]

    scores = []

    for aPlanner in [0, 1]:
        print(f"\n==== Evaluating Planner {aPlanner} ====")
        for i, data in enumerate(problems):
            inputMap, startPos, goalPos = data
            numDOFs = len(startPos.split(","))
            outputSolutionFile = f"tmp_{aPlanner}_{i}.txt"
            startPosString = ",".join(convertPIs(startPos))
            goalPosString = ",".join(convertPIs(goalPos))
            commandPlan = "{} {} {} {} {} {} {}".format(
                executablePath,
                inputMap, numDOFs, startPosString, goalPosString,
                aPlanner, outputSolutionFile)
            commandVerify = "./verifier.out {} {} {} {} {}".format(
                inputMap, numDOFs, startPosString, goalPosString,
                outputSolutionFile)

            try:
                print(f"[Planner {aPlanner}] Problem {i}: Running planner...")
                start = timer()
                subprocess.run(commandPlan.split(" "), check=True)
                timespent = timer() - start
                print(f"[Planner {aPlanner}] Planning done in {timespent:.2f}s")

                print(f"[Planner {aPlanner}] Verifying solution...")
                returncode = subprocess.run(commandVerify.split(" "), check=False).returncode
                if returncode != 0:
                    print("Returned an invalid solution")

                with open(outputSolutionFile) as f:
                    f.readline()
                    solution = [line.split(",")[:-1] for line in f]
                    solution = np.asarray(solution).astype(float)
                    numSteps = solution.shape[0]
                    difsPos = np.abs(solution[1:,] - solution[:-1,])
                    cost = np.minimum(difsPos, np.abs(2*np.pi - difsPos)).sum()
                    success = returncode == 0
                    scores.append([aPlanner, inputMap, i, numSteps, cost, timespent, success])
                    print(f"[Planner {aPlanner}] Steps: {numSteps}, Cost: {cost:.2f}, Success: {success}")

                if success:
                    print(f"[Planner {aPlanner}] Visualizing solution...")
                    commandViz = f"python3 visualizer.py {outputSolutionFile} --gifFilepath=grades_{aPlanner}_{i}.gif --incPrev=1"
                    subprocess.run(commandViz.split(" "), check=True)
                    print(f"[Planner {aPlanner}] Visualization done.")
            except Exception as exc:
                print(f"[Planner {aPlanner}] Failed: {exc} !!")
                scores.append([aPlanner, inputMap, i, -1, -1, 0, False])

    print("\n=== Grading Complete ===")
    df = pd.DataFrame(scores, columns=["planner", "mapName", "problemIndex", "numSteps", "cost", "timespent", "success"])
    df.to_csv(gradingCSV, index=False)
    print(f"Results written to: {gradingCSV}")

            

if __name__ == "__main__":
    graderMain("./planner.out", "test.csv")