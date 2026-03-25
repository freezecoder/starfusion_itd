"""CLI entry point for FusionSQL."""

import click
import yaml
from pathlib import Path

from fusql.utils.logging import setup_logging
from fusql.config import load_config, FusionSQLConfig


@click.group()
@click.option("--log-level", default=None, help="Logging level (overrides config)")
@click.option(
    "--config",
    type=click.Path(exists=False),
    default=None,
    help="Path to config file (default: ./fusql.yaml or ~/.fusql.yaml)",
)
@click.pass_context
def main(ctx: click.Context, log_level: str, config: str):
    """FusionSQL - Gene Fusion Analysis to Microsoft SQL Server."""
    # Load config first
    config_path = Path(config) if config else None
    cfg = load_config(config_path)

    # Store config in context for subcommands
    ctx.obj = {"config": cfg}

    # CLI log-level overrides config
    effective_log_level = log_level if log_level else cfg.log_level
    setup_logging(effective_log_level)


@main.command()
@click.argument("input_dir", type=click.Path(exists=True))
@click.option("--sample-id", help="Override sample ID")
@click.option("--output", "-o", type=click.Path(), help="Output directory")
@click.option("--test-mode", is_flag=True, help="Use CSV output instead of MSSQL")
@click.pass_context
def run(
    ctx: click.Context,
    input_dir: str,
    sample_id: str,
    output: str,
    test_mode: bool,
):
    """Run the full fusion analysis workflow."""
    cfg: FusionSQLConfig = ctx.obj["config"]
    click.echo(f"Running workflow on {input_dir}")
    click.echo(f"Test mode: {test_mode}")
    click.echo(f"Config: MSSQL table={cfg.table_ariba}")


@main.command()
@click.argument("file_path", type=click.Path(exists=True))
@click.option("--output", "-o", type=click.Path(), help="Output file")
@click.pass_context
def parse_ariba(ctx: click.Context, file_path: str, output: str):
    """Parse an Ariba fusion TSV file."""
    from fusql.parsers.ariba import AribaParser
    from pathlib import Path as P

    parser = AribaParser()
    fusions = parser.parse(P(file_path))
    click.echo(f"Parsed {len(fusions)} fusions from {file_path}")

    if output:
        import pandas as pd

        df = pd.DataFrame([f.model_dump() for f in fusions])
        df.to_csv(output, index=False)
        click.echo(f"Written to {output}")


@main.command()
@click.pass_context
def show_config(ctx: click.Context):
    """Display the current configuration."""
    cfg: FusionSQLConfig = ctx.obj["config"]
    click.echo(yaml.dump(cfg.to_dict(), default_flow_style=False, sort_keys=False))


if __name__ == "__main__":
    main()
