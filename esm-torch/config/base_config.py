class BaseConfig:
    # Data paths
    DATA_DIR = 'data/'
    EMBEDDINGS_DIR = 'data/embeddings/'
    RESULTS_DIR = 'experiments/baseline_esm/results/'
    
    # Model parameters
    ESM_DIM = 1280
    HIDDEN_DIMS = [512, 256, 128]
    DROPOUT = 0.3
    
    # Training parameters
    BATCH_SIZE = 32
    LEARNING_RATE = 0.001
    NUM_EPOCHS = 100
    RANDOM_SEED = 42
    
    # Evaluation
    TEST_SIZE = 0.2
    VAL_SIZE = 0.1