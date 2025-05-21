using Documenter
using ParityGate

makedocs(
    sitename = "ParityGate",
    format = Documenter.HTML(),
    modules = [ParityGate],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md"
    ]
)

deploydocs(
    repo = "github.com/yourusername/parity_gate.git"
) 