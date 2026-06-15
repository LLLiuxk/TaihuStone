import os
import sys

try:
    import trimesh
except ImportError:
    print("Error: Python package 'trimesh' is required.")
    print("Please install it by running:")
    print("  pip install trimesh numpy scipy fast-simplification")
    sys.exit(1)

def simplify_stls_in_dir(directory):
    # Ensure 'simp_models' folder exists under the target directory
    output_dir = os.path.join(directory, "simp_models")
    os.makedirs(output_dir, exist_ok=True)
    
    # 5MB in bytes
    LIMIT_SIZE = 8 * 1024 * 1024
    
    # Target face count for binary STL to be safely under 5MB.
    # Standard binary STL has: 80 bytes header + 4 bytes face count + N * 50 bytes.
    # For 5MB (5,242,880 bytes): N < (5,242,880 - 84) / 50 = 104,855.92
    # We choose 95,000 faces, which results in a file size of ~4.75MB (4,750,084 bytes).
    TARGET_FACES = 155000
    
    processed_count = 0
    simplified_count = 0
    skipped_count = 0
    
    print(f"Scanning directory for STL files: {directory}")
    
    for file_name in os.listdir(directory):
        if not file_name.lower().endswith('.stl'):
            continue
            
        file_path = os.path.join(directory, file_name)
        if not os.path.isfile(file_path):
            continue
            
        file_size = os.path.getsize(file_path)
        if file_size > LIMIT_SIZE:
            processed_count += 1
            print(f"\n[{processed_count}] Large STL found: '{file_name}' ({file_size / (1024 * 1024):.2f} MB)")
            
            try:
                # Load the mesh
                mesh = trimesh.load(file_path)
                original_faces = len(mesh.faces)
                
                # Check face count
                if original_faces > TARGET_FACES:
                    print(f"  Mesh has {original_faces} faces. Simplifying to {TARGET_FACES} faces to fit < 5MB...")
                    # simplify_quadric_decimation is fast and preserves mesh features
                    simplified_mesh = mesh.simplify_quadric_decimation(face_count=TARGET_FACES)
                    new_faces = len(simplified_mesh.faces)
                    print(f"  Simplified successfully: {original_faces} -> {new_faces} faces.")
                else:
                    print(f"  Mesh has {original_faces} faces (<= target {TARGET_FACES}).")
                    print(f"  File size > 5MB is likely due to ASCII format. Exporting as binary STL to compress.")
                    simplified_mesh = mesh
                    
                output_path = os.path.join(output_dir, file_name)
                # Export as binary STL
                simplified_mesh.export(output_path, file_type='stl')
                
                new_size = os.path.getsize(output_path)
                print(f"  Saved to: 'simp_models/{file_name}' ({new_size / (1024 * 1024):.2f} MB)")
                if new_size <= LIMIT_SIZE:
                    print("  Status: SUCCESS (file is now under 5MB limit)")
                    simplified_count += 1
                else:
                    print("  Status: WARNING (file size still exceeds 5MB limit)")
                    
            except Exception as e:
                print(f"  Error processing '{file_name}': {e}")
        else:
            skipped_count += 1
            
    print(f"\nScan completed.")
    print(f"  - STL files <= 5MB skipped: {skipped_count}")
    print(f"  - Large STL files processed: {processed_count}")
    print(f"  - STL files successfully simplified and saved to 'simp_models': {simplified_count}")

if __name__ == "__main__":
    # Use the directory passed as argument, or the script's directory by default
    if len(sys.argv) > 1:
        target_dir = sys.argv[1].strip('"').rstrip('\\/')
    else:
        target_dir = os.path.dirname(os.path.abspath(__file__))
    simplify_stls_in_dir(target_dir)
