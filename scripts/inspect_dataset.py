from datasets import load_dataset

try:
    dataset = load_dataset("iidai/codenet", split="train", streaming=True)
    print("Loaded iidai/codenet")
    for sample in dataset:
        print("Keys:", sample.keys())
        print("Sample:", sample)
        break
except Exception as e:
    print(e)
