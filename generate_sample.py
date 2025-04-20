import numpy as np
import random
import subprocess

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

def generate_random_angles(numDOFs=3):
    """ Generates random joint angles between 0 and 2π. """
    # Avoid extreme 0 and 2π
    start = ",".join([f"{random.uniform(0.15, 2.85)}" for _ in range(numDOFs)]) 
    goal = ",".join([f"{random.uniform(0.15, 2.85)}" for _ in range(numDOFs)])
    return convertPIs(start), convertPIs(goal)

def is_valid_configuration(start, goal, map_file, numDOFs):
    """ Runs the C++ verifier to check if the configuration is valid. """
    outputFile = "tmp_check.txt"
    
    # Generate a dummy solution file to pass to the verifier
    with open(outputFile, "w") as f:
        f.write(map_file + "\n")  # First line should be map file
        f.write(start + "\n")     # Start configuration
        f.write(goal + "\n")      # Goal configuration

    # Run the verifier
    commandVerify = f"./verifier.out {map_file} {numDOFs} {start} {goal} {outputFile}"
    try:
        returncode = subprocess.run(commandVerify.split(" "), check=False).returncode
        return returncode == 0  # Return True if configuration is valid
    except:
        return False  # If any error occurs, assume invalid

# Generate and save 5 valid samples
numDOFs = 10  # large_map DOF
numSamples = 5
mapFile = "./map2.txt"
valid_samples = []

while len(valid_samples) < numSamples:
    startPos, goalPos = generate_random_angles(numDOFs)
    startStr, goalStr = ",".join(startPos), ",".join(goalPos)

    if is_valid_configuration(startStr, goalStr, mapFile, numDOFs):
        valid_samples.append(startStr + " " + goalStr)

# Save to fixed_samples.txt
with open("fixed_samples.txt", "w") as f:
    for sample in valid_samples:
        f.write(sample + "\n")

print("Generated 5 VALID start/goal pairs and saved them to fixed_samples.txt.")
