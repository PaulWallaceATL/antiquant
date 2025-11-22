import os
import pandas as pd
from datasets import load_dataset

def ingest_codenet(output_path="data/codenet_sample.csv", num_samples=1000):
    print("Loading dataset...")
    try:
        dataset = load_dataset("iidai/codenet", split="train", streaming=True)
    except Exception as e:
        print(f"Failed to load iidai/codenet: {e}")
        return

    print("Filtering for Python...")
    data = []
    count = 0
    
    for sample in dataset:
        language = sample.get('language', '')
        
        if language == 'Python':
            record = {
                'id': sample.get('id', f'sample_{count}'),
                'language': language,
                'problem_id': sample.get('problem_id', 'unknown'),
                'code': sample.get('code', '')
            }
            
            # Ensure code is not empty
            if record['code']:
                data.append(record)
                count += 1
            
            if count >= num_samples:
                break
    
    print(f"Collected {len(data)} samples.")
    df = pd.DataFrame(data)
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    df.to_csv(output_path, index=False)
    print(f"Saved to {output_path}")

if __name__ == "__main__":
    ingest_codenet()
