"""CLI entry point for FusionSQL."""

import click

from fusql.utils.logging import setup_logging


@click.group()
@click.option("--log-level", default="INFO", help="Logging level")
def main(log_level: str):
    """FusionSQL - Gene Fusion Analysis to Microsoft SQL Server."""
    setup_logging(log_level)


@main.command()
@click.argument("input_dir", type=click.Path(exists=True))
@click.option("--sample-id", help="Override sample ID")
@click.option("--output", "-o", type=click.Path(), help="Output directory")
@click.option("--test-mode", is_flag=True, help="Use CSV output instead of MSSQL")
def run(input_dir: str, sample_id: str, output: str, test_mode: bool):
    """Run the full fusion analysis workflow."""
    click.echo(f"Running workflow on {input_dir}")
    click.echo(f"Test mode: {test_mode}")


@main.command()
@click.argument("file_path", type=click.Path(exists=True))
@click.option("--output", "-o", type=click.Path(), help="Output file")
def parse_ariba(file_path: str, output: str):
    """Parse an Ariba fusion TSV file."""
    from fusql.parsers.ariba import AribaParser
    from pathlib import Path
    
    parser = AribaParser()
    fusions = parser.parse(Path(file_path))
    click.echo(f"Parsed {len(fusions)} fusions from {file_path}")
    
    if output:
        import pandas as pd
        df = pd.DataFrame([f.model_dump() for f in fusions])
        df.to_csv(output, index=False)
        click.echo(f"Written to {output}")


if __name__ == "__main__":
    main()
