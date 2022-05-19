from torch.utils.data import Dataset

class ShortReadsDataset(Dataset):
    def __init__(self, fasta_file, ) -> None:
        super().__init__()