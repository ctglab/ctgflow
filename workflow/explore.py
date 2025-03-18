from pathlib import Path

from snakemake.api import (
    OutputSettings,
    ResourceSettings,
    SnakemakeApi,
    StorageSettings,
)

with SnakemakeApi(
    OutputSettings(
        verbose=False,
        show_failed_logs=True,
    ),
) as snakemake_api:
    workflow_api = snakemake_api.workflow(
        storage_settings=StorageSettings(),
        resource_settings=ResourceSettings(),
        snakefile=Path("Snakefile"),
    )
    dag_api = workflow_api.dag()
    # Go on by calling methods of the dag api.
