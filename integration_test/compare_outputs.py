def compare_files(file1, file2):
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        for line1, line2 in zip(f1, f2):
            if line1 != line2:
                print(f"Difference found:\n{file1}: {line1}\n{file2}: {line2}")
                return False
    return True

if __name__ == "__main__":
    if compare_files('positions1.xyz', 'positions2.xyz'):
        print("Test Passed: Files are identical.")
    else:
        print("Test Failed: Files differ.")