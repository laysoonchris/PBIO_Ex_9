# s31774_2026.py
# Numer albumu: s31774
# Data: 2026-05-07
# Opis: Generator losowych sekwencji DNA z zapisem w formacie FASTA.
# Obsługuje: tryb wsadowy, konfigurowalny rozkład nukleotydów,
# transkrypcję in silico (T→U) oraz analizę okien przesuwnych GC
# z generowaniem wykresu liniowego.

import random
import csv
import re
from typing import Optional

# matplotlib jest potrzebny tylko do wykresu dla podpunktu 8
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

def validate_positive_int(prompt: str, min_val: int = 1, max_val: int = 100_000) -> int:
    """Pobiera od użytkownika liczbę całkowitą z zakresu.
    W przypadku błędu powtarza pytanie.
    """
    while True:
        raw_input = input(prompt)
        try:
            value = int(raw_input)
            if min_val <= value <= max_val:
                return value
            print(f"Błąd: wartość musi być liczbą całkowitą z zakresu [{min_val}, {max_val}].")
        except ValueError:
            print(f"Błąd: wartość musi być liczbą całkowitą z zakresu [{min_val}, {max_val}].")

def validate_seq_id(prompt: str) -> str:
    """Pobiera identyfikator sekwencji zgodny ze specyfikacją FASTA.
    ID nie może być puste ani zawierać białych znaków. Pyta ponownie przy błędzie.
    """
    while True:
        seq_id = input(prompt).strip()
        if seq_id and not re.search(r"\s", seq_id):
            return seq_id
        print("Błąd: ID nie może być puste ani zawierać białych znaków.")

def generate_sequence(length: int) -> str:
    """Zwraca losową sekwencję DNA o zadanej długości."""
    return "".join(random.choices(['A', 'C', 'G', 'T'], k=length))

def _generate_sequence_weighted(length: int, weights: dict) -> str:
    """Zwraca losową sekwencję DNA z niestandardowym rozkładem nukleotydów."""
    nucleotides = ['A', 'C', 'G', 'T']
    w = [weights["A"], weights["C"], weights["G"], weights["T"]]
    return "".join(random.choices(nucleotides, weights=w, k=length))

def format_fasta(seq_id: str, description: str, sequence: str, line_width: int = 80) -> str:
    """Zwraca sformatowany rekord FASTA jako string."""
    header = f">{seq_id} {description}" if description else f">{seq_id}"
    lines = [sequence[i: i + line_width] for i in range(0, len(sequence), line_width)]
    return header + "\n" + "\n".join(lines) + "\n"

def calculate_stats(sequence: str) -> dict:
    """Zwraca słownik ze statystykami sekwencji.
    Klucze: "A", "C", "G", "T" (wartości float, %),
            "GC" (wartość float, %).
    """
    bio_seq = "".join(c for c in sequence if c in "ACGT")
    total_length = len(bio_seq)

    if total_length == 0:
        return {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0, "GC": 0.0}

    counts = {nucleotide: bio_seq.count(nucleotide) for nucleotide in "ACGT"}
    stats = {nucleotide: round(counts[nucleotide] / total_length * 100, 2) for nucleotide in "ACGT"}
    stats["GC"] = round((counts['G'] + counts['C']) / total_length * 100, 2)
    return stats


def print_stats(stats: dict, length: int) -> None:
    """Wypisuje statystyki sekwencji w czytelnym formacie na stdout."""
    print(f"\nStatystyki sekwencji (n={length}):")
    for nucleotide in ["A", "C", "G", "T"]:
        print(f"  {nucleotide}: {stats[nucleotide]:.2f}%")
    print(f"  GC-content: {stats['GC']:.2f}%")


def insert_name(sequence: str, name: str) -> str:
    """Wstawia imię w losową pozycję sekwencji.
    Imię zapisane małymi literami.
    """
    position = random.randint(0, len(sequence))
    return sequence[:position] + name.lower() + sequence[position:]

# feature 1: batch mode

def generate_batch(count: int, length: int, base_id: str, description: str, name: str, weights: Optional[dict], output_file: str) -> None:
    """Generuje wiele sekwencji DNA i zapisuje je do jednego pliku multi-FASTA.
    ID kolejnych sekwencji: {base_id}_{numer:03d} (np. Seq_001, Seq_002).
    Statystyki liczone na czystej sekwencji (bez imienia).
    """
    with open(output_file, "w", encoding="utf-8") as f:
        for i in range(1, count + 1):
            seq_id = f"{base_id}_{i:03d}"
            seq = _generate_sequence_weighted(length, weights) if weights else generate_sequence(length)
            seq_with_name = insert_name(seq, name)
            
            record = format_fasta(seq_id, description, seq_with_name)
            f.write(record)

            stats = calculate_stats(seq)
            print(f"[{i}/{count}] {seq_id}  GC={stats['GC']:.2f}%")

    print(f"\nZapisano {count} sekwencji do pliku: {output_file}")

# feature 2: konfigurowalny rozkład nukleotydów

def get_nucleotide_distribution() -> dict:
    """Pobiera od użytkownika procentowy udział każdego nukleotydu.
    Waliduje że suma wynosi 100 (tolerancja 0.01) i wszystkie wartości są nieujemne.
    Zwraca słownik {'A': float, 'C': float, 'G': float, 'T': float}.
    """
    print("Podaj procentowy udział nukleotydów (suma musi wynosić 100):")
    while True:
        try:
            a = float(input("  A (%): "))
            c = float(input("  C (%): "))
            g = float(input("  G (%): "))
            t = float(input("  T (%): "))
            total = a + c + g + t

            if abs(total - 100.0) < 0.01 and all(v >= 0 for v in [a, c, g, t]):
                return {"A": a, "C": c, "G": g, "T": t}

            print(f"Błąd: suma wynosi {total:.2f}, musi wynosić 100. Spróbuj ponownie.")
        except ValueError:
            print("Błąd: podaj liczby dziesiętne (np. 30.5).")

# feature 5: transkrypcja in silico

def transcribe(dna_sequence: str) -> str:
    """Generuje sekwencję mRNA z sekwencji DNA przez zamianę T -> U.
    Działa na czystej sekwencji bez imienia — małe litery pozostają niezmienione.
    """
    return dna_sequence.replace("T", "U")

def sliding_window_gc(sequence: str, window_size: int) -> list:
    """Oblicza GC-content w oknie przesuwnym po sekwencji.
    Zwraca listę słowników z kluczami 'start_position' (1-based) i 'gc_content' (%).
    Małe litery (imię) są filtrowane przed obliczeniami.
    """
    bio_seq = "".join(c for c in sequence if c in "ACGT")
    gc_content_values = []

    for i in range(len(bio_seq) - window_size + 1):
        window = bio_seq[i: i + window_size]
        gc = (window.count("G") + window.count("C")) / window_size * 100
        gc_content_values.append({"start_position": i + 1, "gc_content": round(gc, 2)})
        
    return gc_content_values


def save_gc_content_to_csv(gc_content_values: list, filename: str) -> None:
    """Zapisuje wyniki analizy okien przesuwnych do pliku CSV.
    Kolumny: start_position, gc_content.
    """
    with open(filename, mode='w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=["start_position", "gc_content"])
        writer.writeheader()
        writer.writerows(gc_content_values)
        
    print(f"Wyniki sliding window zapisane do: {filename}")

# feature 8: GC-content chart

def plot_gc_chart(results: list, window_size: int, seq_id: str) -> str:
    """Generuje wykres liniowy GC-content ze sliding window i zapisuje do PNG.
    Linia przerywana oznacza średni GC-content całej sekwencji.
    Zwraca ścieżkę do pliku lub pusty string gdy matplotlib niedostępny.
    """
    if not MATPLOTLIB_AVAILABLE:
        print("Uwaga: matplotlib niedostępny.")
        
        return ""

    positions = [r["start_position"] for r in results]
    gc_values = [r["gc_content"] for r in results]

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(positions, gc_values, linewidth=0.8, color="#2196F3")

    mean_gc = sum(gc_values) / len(gc_values)
    ax.axhline(mean_gc, color="#F44336", linestyle="--", linewidth=1, label=f"Średni GC: {mean_gc:.1f}%")

    ax.set_xlabel("Pozycja startu okna (nt)")
    ax.set_ylabel("GC-content (%)")
    ax.set_title(f"GC-content – okno przesuwne (szerokość={window_size} nt) | {seq_id}")
    ax.set_ylim(0, 100)
    ax.legend()
    fig.tight_layout()

    output_path = f"{seq_id}_gc_chart.png"
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"Wykres GC-content zapisany do: {output_path}")
    return output_path

def main():
    """Wiadomo."""
    print("Generator sekwencji DNA (format FASTA)\n")

# Krok 1: zapytaj o tryby działania
    batch_mode = input("Tryb wsadowy (batch)? (t/n): ").strip().lower() == "t"
    custom_dist = input("Konfigurowalny rozkład nukleotydów? (t/n): ").strip().lower() == "t"

    weights = get_nucleotide_distribution() if custom_dist else None

# Krok 2: parametry wspólne dla obu trybów
    length = validate_positive_int("Podaj długość sekwencji: ")
    seq_id = validate_seq_id("Podaj ID sekwencji: ")
    description = input("Podaj opis sekwencji: ")
    name = input("Podaj imię: ").strip()

# TRYB WSADOWY
    if batch_mode:
        count = validate_positive_int("Podaj liczbę sekwencji: ", min_val=1, max_val=1000)
        output_file = f"{seq_id}_batch.fasta"
        generate_batch(count, length, seq_id, description, name, weights, output_file)

# TRYB POJEDYNCZY
    else:
        # Generowanie sekwencji
        sequence = _generate_sequence_weighted(length, weights) if weights else generate_sequence(length)
        # Oblicz statystyki zanim wstawimy imię
        stats = calculate_stats(sequence)
        # Wstaw imię małymi literami w losową pozycję
        sequence_with_name = insert_name(sequence, name)

        # feature 5
        mrna_sequence = transcribe(sequence)
        mrna_id = f"{seq_id}_mRNA"

        # save to FASTA file
        output_file = f"{seq_id}.fasta"
        with open(output_file, "w", encoding="utf-8") as f:
            # Rekord 1: oryginalna sekwencja DNA z imieniem
            f.write(format_fasta(seq_id, description, sequence_with_name))
            # Rekord 2: sekwencja mRNA (transkrypt)
            f.write(format_fasta(mrna_id, f"mRNA transcript | {description}", mrna_sequence))

        print(f"\nSekwencja zapisana do pliku: {output_file}")
        print_stats(stats, length)

        # feature 8 analiza okien przesuwnych i wykres GC
        sw_choice = input("\nWykres GC-content (sliding window)? (t/n): ").strip().lower()
        if sw_choice == "t":
            max_window = min(length, 100_000)
            window_size = validate_positive_int(
                "Podaj szerokość okna (nt): ", min_val=1, max_val=max_window
            )

            sw_results = sliding_window_gc(sequence, window_size)

            # zapisz csv
            csv_file = f"{seq_id}_gc_window.csv"
            save_gc_content_to_csv(sw_results, csv_file)

            # generowanie wykresu PNG
            plot_gc_chart(sw_results, window_size, seq_id)

if __name__ == "__main__":
    main()